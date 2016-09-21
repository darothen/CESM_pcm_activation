//
// Based on: http://jsfiddle.net/QDP9H/
//


var json_file = "mode.json";
var modes_array = [
  "ACC", "AIT", "MOS", "MBS", "DST01", "DST02", "SSLT01"
];

// Breadcrumb dimensions: width, height, spacing, width of tip/tail.
var b = {
    w: 75, h: 30, s: 3, t: 10
};

var width = 600,
    height = 600,
    offset = 50,
    radius = (Math.min(width, height) / 2) - offset,
    color = d3.scale.ordinal()
      .domain(modes_array)
      .range(colorbrewer.Dark2[modes_array.length]),
    totalCount = 0.0,
    notchRads = Math.PI/8;
    iterLevels = 0;

var svg = d3.select("#chart").append("svg:svg")
    .attr("width", width)
    .attr("height", height)
    .append("svg:g")
    .attr("id", "container")
    .attr("transform", "translate(" + width / 2 + "," + height * .52 + ")");

var partition = d3.layout.partition()
    .sort(null)
    .size([2 * Math.PI - notchRads, 100])
    .value(function(d) {
      totalCount += d.size; // update total number of samples
      if ( iterLevels < d.depth ) {
        iterLevels += 1;
      }
      return d.size;
    });

var arc = d3.svg.arc()
    .startAngle(function(d) { return d.x + notchRads/2; })
    .endAngle(function(d) { return d.x + d.dx + notchRads/2; })
    .innerRadius(function(d) { return radius * (d.y) / 100; })
    .outerRadius(function(d) { return radius * (d.y + d.dy) / 100; });

var legend = null;

// This main method creates the visualization as a callback
// to the JSON-reading method
d3.json(json_file, function(error, root) {
  if (error) throw error;

  console.log(root);

  // Setup breadcumb trail elements
  initializeBreadcrumbTrail();

  var path = svg.datum(root).selectAll("path")
      .data(partition.nodes)
      .enter().append("path")
        .attr("class", "arc")
        .attr("name", function(d) { return d.name })
        .attr("display", function(d) { return d.depth ? null : "none"; }) // hide inner ring
        .attr("d", arc)
        .style("fill", function(d) {
          return color(d.name);
        })
        .on("mouseover", function(d) {
            console.log(d.name, d.depth, d.x, d.dx);
            mouseover(d);
        });

    // Add the mouseleave handler to the bounding circle
    d3.select("#container").on("mouseleave", mouseleave);

    // Add text inside the 'notch' denoting the iterations);
    var g = svg.selectAll("g")
      .data(range(1, iterLevels))
      .enter().append("svg:g");

    g.append("svg:text")
      .attr("y", function(d) {
        var vert = height/2;
        var levelWidth = radius/(iterLevels + 1.0);
        return -vert + levelWidth * ((iterLevels - d) + 0.5) + offset;
      })
      .attr("dy", "1em")
      .attr("text-anchor", "middle")
      .style("font-weight", 800)
      .text(function(d) { return d });


    svg.insert("svg:g", ":first-child")
      .attr("transform", function(d) {
          var levelWidth = radius/(iterLevels + 1.0);
          var titleLoc = -height/2 - levelWidth/10. + offset
          return "translate(0," + titleLoc + ")";
      })
      .append("svg:text")
        .attr("class", "title")
        .attr("dy", "0.35em")
        .attr("text-anchor", "middle")
        .style({
         "fill": "#000",
         "font-weight": "bold"
        })
        .text("Iteration");

    drawLegend();
});

// Fade all but the current sequence, and show it in the
// breadcrumb trail
function mouseover(d) {
  var percentage = (100 * d.value / totalCount).toPrecision(3);
  var percentageString = percentage + "%";
  if ( percentage < 0.1 ) {
    percentageString = "< 0.1%";
  }

  d3.select("#percentage")
    .text(percentageString);

  var sequenceArray = getAncestors(d);
  updateBreadcrumbs(sequenceArray, percentageString);

  // Fade all the segments
  d3.selectAll("path")     // arcs in sunburst
    .style("opacity", 0.3);
  d3.selectAll(".legend_rec")     // labels in legend
    .style("opacity", 0.3);

  // Highlight only those that are an ancestor of
  // the current segment
  svg.selectAll("path")
    .filter(function(node) {
      return (sequenceArray.indexOf(node) >= 0);
    })
    .style("opacity", 1);
}

// Restore everything to full opacity when moving off the vis
function mouseleave(d) {

  // Hide the breadcrumb trail
  d3.select("#trail")
    .style("visibility", "hidden");

  // Deactivate all segments during transition
  d3.selectAll("path").on("mouseover", null);

  // Transition each segment to full opacity then reactivate it
  d3.selectAll("path")
    .transition()
    .duration(0)
    .style("opacity", 1)
    .each("end", function() {
      d3.select(this).on("mouseover", mouseover);
    });
  d3.selectAll(".legend_rec")     // labels in legend
    .transition()
    .duration(0)
    .style("opacity", 1);
}

// Function for drawing the legend
function drawLegend() {

  // Dimensions of legend item: width, height, spacing, radius of rounded rect.
  var li = {
    w: 75, h: 30, s: 3, r: 3
  };

  legend = d3.select("#legend").append("svg:svg")
    .attr("width", li.w)
    .attr("height", color.domain().length * (li.h + li.s));

  var g = legend.selectAll("g")
    .data(modes_array)
    .enter().append("svg:g")
    .attr("transform", function(d, i) {
        return "translate(0," + (i+1) * (li.h + li.s) + ")";
    });

  g.append("svg:rect")
    .attr("class", "legend_rec")
    .attr("mode", function(d) { return d; })
    .attr("rx", li.r)
    .attr("ry", li.r)
    .attr("width", li.w)
    .attr("height", li.h)
    .style("fill", function(d) { return color(d); });

  g.append("svg:text")
    .attr("x", li.w / 2)
    .attr("y", li.h / 2)
    .attr("dy", "0.35em")
    .attr("text-anchor", "middle")
    .text(function(d) { return d });

  legend.insert("svg:g", ":first-child")
    .attr("transform", function(d) {
        return "translate(0,0)";
    })
    .append("svg:text")
      .attr("class", "title")
      .attr("x", li.w / 2)
      .attr("y", li.h / 2)
      .attr("dy", "0.35em")
      .attr("text-anchor", "middle")
      .style({
       "fill": "#000",
       "font-weight": "bold"
      })
      .text("Mode Name");
}

// Breadcrumb / sequence analysis
function initializeBreadcrumbTrail() {
  // Add the svg area for the graphic
  var trail = d3.select("#sequence").append("svg:svg")
    .attr("width", width)
    .attr("height", 50)
    .attr("id", "trail");
  // Add the label at the end, for the percentage
  trail.append("svg:text")
    .attr("id", "endlabel")
    .style("fill", "#000");

}

// Generate a string that describes the points of a breadcrumb
// polygon
function breadcrumbPoints(d, i) {
  var points = [];
  points.push("0,0"); // these vertices are for notched rectangles
  points.push(b.w + ",0");
  points.push(b.w + b.t + "," + (b.h / 2));
  points.push(b.w + "," + b.h);
  points.push("0," + b.h);
  if (i > 0) { // Leftmost breadcrumb; don't include 6th vertex.
    points.push(b.t + "," + (b.h / 2));
  }
  return points.join(" ");
}

// Update the breadcumb trail sequence and percentage
function updateBreadcrumbs(nodeArray, percentageString) {

  // Data join; key function combines name and depth
  // (= position in sequence)
  var g = d3.select("#trail").selectAll("g")
    .data(nodeArray, function(d) { return d.name + d.depth; });

  // Add breadcrumb and label for entering nodes
  var entering = g.enter().append("svg:g");
  entering.append("svg:polygon")
    .attr("points", breadcrumbPoints)
    .style("fill", function(d) { return color(d.name); });
  entering.append("svg:text")
    .attr("x", (b.w + b.t) / 2)
    .attr("y", b.h/2)
    .attr("dy", "0.35em")
    .attr("text-anchor", "middle")
    .text(function(d) { return d.name; });

  // Set position for entering and updating nodes
  g.attr("transform", function(d, i) {
    return "translate(" + i*(b.w + b.s) + ",0)";
  });

  // Remove exiting nodes.
  g.exit().remove()

  // Now move and update the percentage at the end.
  d3.select("#trail").select("#endlabel")
    .attr("x", (nodeArray.length + 0.5) * (b.w + b.s))
    .attr("y", b.h / 2)
    .attr("dy", "0.35em")
    .attr("text-anchor", "middle")
    .text(percentageString);

  // Make the breadcrumb trail visible if it's hidden
  d3.select("#trail")
    .style("visibility", "");
}

// Given a node in a partition layout, return an array of all of
// its ancestor nodes, highest first, but excluding the root.
function getAncestors(node) {
  var path = [];
  var current = node;
  while ( current.parent ) {
    path.unshift(current);
    current = current.parent;
  }
  return path;
}

// Range function
function range(start, end) {
  var rangeArray = [];
  for (var i=start; i <= end; i++) {
    rangeArray.push(i);
  }
  return rangeArray;
}

d3.select(self.frameElement).style("height", (height+100) + "px");
