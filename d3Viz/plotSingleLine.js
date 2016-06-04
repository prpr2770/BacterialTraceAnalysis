function plotSingleLine(linePlot_data, linePlot_name, divID, margin, width, height){
/*
// =============================================================================
// CSS Style and the initial code to run the script!
// =============================================================================

var cellData = [
{cellID: 1, voltSig : [1, 2, 3, 4, 5,1, 2, 3, 4, 5], caSig: [1, 2, 3, 4, 5, 1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5], xPos: 1, yPos: 1 },
{cellID: 2, voltSig : [1, 2, 3, 4, 5,5, 4, 3, 2, 1], caSig: [1, 2, 3, 4, 5, 1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5], xPos: 2, yPos: 2 },
{cellID: 3, voltSig : [1, 2, 3, 4, 5,6, 7, 8, 7, 6], caSig: [1, 2, 3, 4, 5, 1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5], xPos: 3, yPos: 3 },
{cellID: 4, voltSig : [1, 2, 3, 4, 5,5, 6, 7, 8, 9], caSig: [1, 2, 3, 4, 5, 1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5], xPos: 4, yPos: 4 },
{cellID: 5, voltSig : [1, 2, 3, 4, 5,4, 3, 2, 1, 0], caSig: [1, 2, 3, 4, 5, 1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5], xPos: 5, yPos: 5 }];

var margin = {top: 20, right: 80, bottom: 30, left: 50},
width = 960 - margin.left - margin.right,
height = 500 - margin.top - margin.bottom;

var divID = "#multilineplot";

var linePlot_data = cellData[0].voltSig;
plotSingleLine(linePlot_data, divID, margin, width, height);


<style>
  path {
    stroke: steelblue;
    stroke-width: 1;
    fill: none;
  }

  .axis {
    shape-rendering: crispEdges;
  }

  .x.axis line {
    stroke: lightgrey;
  }

  .x.axis .minor {
    stroke-opacity: .5;
  }

  .x.axis path {
    display: none;
  }

  .y.axis line, .y.axis path {
    fill: none;
    stroke: #000;
  }

  .selected {
    fill: blue;
  }

  .brush .extent {
  fill-opacity: .1;
  stroke: #fff;
  shape-rendering: crispEdges;
}

</style>



// =============================================================================
*/

var xScale = d3.scale.linear().range([0, width]);
var yScale = d3.scale.linear().range([height, 0]);


  // Add an SVG element with the desired dimensions and margin.
  var graph_width = width + margin.right + margin.left;
  var graph_height = height + margin.top + margin.bottom;
  var graph = d3.select(divID).append("svg:svg")
        .attr("class","aGraph")
        .attr("width", graph_width)
        .attr("height", graph_height)
        .classed("singleline",true)
        .append("svg:g") // what does this do?
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
        .attr("class","aLinePlot");

graph.data(function(){
  xScale.domain([0, linePlot_data.length]);
  yScale.domain([d3.min(linePlot_data), d3.max(linePlot_data)]);
  return [linePlot_data];
});


// create a line function that can convert data[] into x and y points
var linePlot = d3.svg.line()
  .x(function(datum,indx) {  return xScale(indx); })
  .y(function(datum) { return yScale(datum);});

  // create yAxis
  var xAxis = d3.svg.axis().scale(xScale).tickSize(-height).tickSubdivide(true);
  // Add the x-axis.
  graph.append("svg:g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + height+ ")")
        .call(xAxis);


  // create left yAxis
  var yAxisLeft = d3.svg.axis().scale(yScale).ticks(4).orient("left");
  // Add the y-axis to the left
  graph.append("svg:g")
        .attr("class", "y axis")
        .attr("transform", "translate(-10,0)")
        .call(yAxisLeft);

    // Add the line by appending an svg:path element with the data line we created above
  // do this AFTER the axes above so that the line is above the tick-lines

      graph.append("svg:path")
        .attr("class","view_plot")
        .attr("d", function(d){return linePlot(d);});

        graph.append("text")// title of cell, place it near last datapoint!
                .attr("class", "view_plot_name")
                  .attr("transform", function(d) { var indx_last_element = d.length -1; return "translate(" + xScale(indx_last_element) + "," + yScale(d[indx_last_element]) + ")"; })
                  .attr("x", 3)
                  .attr("dy", ".35em")
                  .text(function() { return linePlot_name; });


return [xScale, yScale, xAxis, yAxisLeft, linePlot];

}
