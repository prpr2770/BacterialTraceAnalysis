data = [{x:1,y:1}, {x:2,y:2}, {x:3,y:3}, {x:4,y:4}];
width = 400;
height = 300;

var xScale = d3.scale.linear().domain([0,5]).range([0,width]);
var yScale = d3.scale.linear().domain([0,5]).range([height,0]);

var divdots = d3.select(".divdots");
var dots = divdots.append("svg:svg")
              .attr("class","dots")
              .attr("width",width)
              .attr("height", height);

// =============================================================================
var brush = d3.svg.brush()
	.x(d3.scale.identity().domain([0, width]))
	.y(d3.scale.identity().domain([0, height]))
	//.on("brushstart", brushstart)
	.on("brush", brushmove)
	.on("brushend", brushend);

function brushmove() {
          var e = brush.extent();
          dots.selectAll(".g_cell").classed("selected", function(d) {
              return !( e[0][0] > xScale(d.x) || xScale(d.x)> e[1][0] || e[0][1] > yScale(d.y) || yScale(d.y) > e[1][1]);
              });
        }

// If the brush is empty, select all circles.
function brushend() {
          if (brush.empty()) { // When NO AREA is selected!
            //alert("Brush: EMPTY");
            // update class attribute of the cells.
            dots.selectAll(".selected").classed("selected", false);
            }
      }

// assign classname "brush" to the g_brush contained in divdots
dots.classed("brush",true).call(brush);

// =============================================================================
// Defining circles


var dot = dots.selectAll(".g_cell").data(data);

dot.enter()
      .append("svg:g")
        .attr("class","g_cell")
        .append("circle")
          .attr("class","cell")
          .attr("cx", function(d) { return xScale(d.x); })
          .attr("cy", function(d) { return yScale(d.y); })
          .attr("r", 0)
          .transition()
            .attr("r",5);
