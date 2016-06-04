function plotMultiLine(myData, divID, margin, width, height) {
/*
// =============================================================================
// CSS Style and the initial code to run the script!
// =============================================================================
var myData = [
{cellID: 1, voltSig : [1, 2, 3, 4, 5,1, 2, 3, 4, 5], caSig: [1, 2, 3, 4, 5, 1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5], xPos: 1, yPos: 1 },
{cellID: 2, voltSig : [1, 2, 3, 4, 5,5, 4, 3, 2, 1], caSig: [1, 2, 3, 4, 5, 1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5], xPos: 2, yPos: 2 },
{cellID: 3, voltSig : [1, 2, 3, 4, 5,6, 7, 8, 7, 6], caSig: [1, 2, 3, 4, 5, 1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5], xPos: 3, yPos: 3 },
{cellID: 4, voltSig : [1, 2, 3, 4, 5,5, 6, 7, 8, 9], caSig: [1, 2, 3, 4, 5, 1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5], xPos: 4, yPos: 4 },
{cellID: 5, voltSig : [1, 2, 3, 4, 5,4, 3, 2, 1, 0], caSig: [1, 2, 3, 4, 5, 1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5], xPos: 5, yPos: 5 }];

var margin = {top: 20, right: 80, bottom: 30, left: 50},
width = 960 - margin.left - margin.right,
height = 500 - margin.top - margin.bottom;

var divID = "#multilineplot";
plotMultiLine(myData, divID, margin, width, height);

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

    .legend rect {
      fill:white;
      stroke:black;
      opacity:0.8;}

  }

</style>
// =============================================================================
*/
var cellIDs = [];
for (var i = 0, len = myData.length; i < len; i++) {
  cellIDs.push(myData[i].cellID);
}

var color = d3.scale.category10().domain(cellIDs);



var xScale = d3.scale.linear().range([0, width]);
var yScale = d3.scale.linear().range([height, 0]);
var svg = d3.select("body").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");


// Add an SVG element with the desired dimensions and margin.
var graph_width = width + margin.left + margin.right;
var graph_height = height + margin.top + margin.bottom;
var graph = d3.select("#"+divID).append("svg:svg")
      .attr("class","aGraph_"+divID)
      .classed(divID,true)
      .attr("width", graph_width)
      .attr("height", graph_height)
      .append("svg:g") // what does this do?
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
      .attr("class","manyLinePlot");


graph.data(function(){
  // update the Scales!
  //console.log(myData[0]);
  xScale.domain([0, myData[0].voltSig.length]);
  yScale.domain([
    d3.min(myData,function(c){return d3.min(c.voltSig);}),
    d3.max(myData,function(c){return d3.max(c.voltSig);})
  ]);

  //--------------------------------------------------------------------------
  // associate Data_set to the graph!
  return myData; // assosciate subset of myData to the plot_graph.
  });

//--------------------------------------------------------------------------
// set the axis
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

  //----------------------------------------------------------------------------
  // plot the voltage signals
var sig = graph.selectAll(".sig")
      .data(myData)
        .enter().append("g")
          .attr("class", "sig");

sig.append("path")
    .attr("class", "plotline")
    .attr("d", function(d) { return linePlot(d.voltSig); })
    .attr("data-legend",function(d) { return "cell_"+d.cellID; })
    .style("stroke", function(d) { return color(d.cellID); });

sig.append("text")// title of cell, place it near last datapoint!
        .attr("class", "plotname")
          .attr("transform", function(d) { var indx_last_element = d.voltSig.length -1; return "translate(" + xScale(indx_last_element) + "," + yScale(d.voltSig[indx_last_element]) + ")"; })
          .attr("x", 3)
          .attr("dy", ".35em")
          .text(function(d) { return "cell_"+d.cellID; });


//==============================================================================
///*
// Create the layer in which legend shall be displayed.
sig.append("g")
            .attr("class","legend")
            .attr("transform","translate("+ (width-margin.left-margin.right) +",0)")
            .style("font-size","10px")
            .call(d3.legend);
  //*/




}
