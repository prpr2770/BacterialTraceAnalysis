function plotMultiLine6(myData, params) {

console.log(myData);
console.log(params);
var divID = params.divID;
var class_name = params.classname;
var margin = params.margin;
var width = params.width;
var height = params.height;

console.log(margin);
console.log(margin.left);
console.log(margin.right);
console.log(margin.top);
console.log(margin.bottom);



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
var graph = d3.select("."+class_name).selectAll("#"+divID).append("svg:svg")
      .attr("class","aGraph_"+divID+"_"+class_name)
      .classed(divID,true)
      .attr("width", graph_width)
      .attr("height", graph_height)
      .append("svg:g") // what does this do?
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
      .attr("class","manyLinePlot");


graph.data(function(){
  // update the Scales!
  //console.log(myData[0]);

  switch(class_name){
    case "volt":
          xScale.domain([0, myData[0].voltSig.length]);
          yScale.domain([
            d3.min(myData,function(c){return d3.min(c.voltSig);}),
            d3.max(myData,function(c){return d3.max(c.voltSig);})
          ]);
          break;
    case "calc":
          xScale.domain([0, myData[0].caSig.length]);
          yScale.domain([
            d3.min(myData,function(c){return d3.min(c.caSig);}),
            d3.max(myData,function(c){return d3.max(c.caSig);})
          ]);
          break;
  }



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


          switch(class_name){
            case "volt":
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
                  break;


            case "calc":
                  sig.append("path")
                      .attr("class", "plotline")
                      .attr("d", function(d) { return linePlot(d.caSig); })
                      .attr("data-legend",function(d) { return "cell_"+d.cellID; })
                      .style("stroke", function(d) { return color(d.cellID); });

                  sig.append("text")// title of cell, place it near last datapoint!
                          .attr("class", "plotname")
                            .attr("transform", function(d) { var indx_last_element = d.caSig.length -1; return "translate(" + xScale(indx_last_element) + "," + yScale(d.caSig[indx_last_element]) + ")"; })
                            .attr("x", 3)
                            .attr("dy", ".35em")
                            .text(function(d) { return "cell_"+d.cellID; });
                  break;
          }




//==============================================================================
// Create the layer in which legend shall be displayed.
sig.append("g")
            .attr("class","legend")
            .attr("transform","translate("+ (width-margin.left-margin.right) +",0)")
            .style("font-size","10px")
            .call(d3.legend);





}
