function plotPetridish(linkData, cellData, params){

// ------------------------------------------------------------------
var edge_type = "arc"; // "arc", "line", "curved"
var max_edge_thickness = 1;

var pd_width = params.pd_width;
var pd_height = params.pd_height;
var pd_margin = params.pd_margin;
var pd_classname = params.pd_class_name;
var pd_divID = params.pd_divID;

var petriWidth = pd_width - pd_margin.left - pd_margin.right;
var petriHeight = pd_height - pd_margin.top - pd_margin.bottom;
var petridish = d3.select("."+pd_classname).select("#"+pd_divID).attr("width", pd_width ).attr("height", pd_height );

///*
// =====================================================================================================
// Defining Scales to be Used for spatial distribution
var pd_xScale = d3.scale.linear().range([pd_margin.left, petriWidth]);
var pd_yScale = d3.scale.linear().range([petriHeight, pd_margin.top]);

// set domain by considering dataset.
pd_xScale.domain([d3.min(cellData, function(d){return d.xPos;}), d3.max(cellData, function(d){return d.xPos;})]);
pd_yScale.domain([d3.min(cellData, function(d){return d.yPos;}), d3.max(cellData, function(d){return d.yPos;})]);


// =====================================================================================================
// Defining and Using Brush
var brush = d3.svg.brush()
.x(d3.scale.identity().domain([0, pd_width]))
.y(d3.scale.identity().domain([0, pd_height]))
//.on("brushstart", brushstart)
.on("brush", brushmove)
.on("brushend", brushend);

// Highlight the selected circles.
function brushmove() {
var e = brush.extent();
petridish.selectAll(".g_cell").classed("selected", function(d) {
return !( e[0][0] > pd_xScale(d.xPos) || pd_xScale(d.xPos) > e[1][0] || e[0][1] > pd_yScale(d.yPos) || pd_yScale(d.yPos) > e[1][1]);
});

}

// If the brush is empty, select all circles.
function brushend() {
      // remove the display of cell_names always!
      petridish.selectAll(".cellname").remove();
      if (brush.empty()) { // When NO AREA is selected!
        //alert("Brush: EMPTY");
        // update class attribute of the cells.
        petridish.selectAll(".selected").classed("selected", false);
}

// -------------------------------------------------------------------
// display cellID of selectedCells
//var nodes_g = petridish.selectAll(".g_cell");
var selected_nodes = petridish.selectAll(".selected");

// if NO CELLS are SELECTED:
if (selected_nodes.node() === null){
//alert("selected_nodes.node() == null");
// remove the names of previously selected cells.
selected_nodes.selectAll(".cellname").remove();
selected_nodes.classed("selected", false);
  }
else{
//alert("selected_nodes.node() != null");

selected_nodes.append("text")// title of cell, place it near center of Cell
.attr("class", "cellname")
//.attr("transform", function(d) { return "translate(" + d.cx+ "," + d.cy+ ")"; })
.attr("x", function(d){return d.cx;})
.attr("y", function(d){return d.cy;})
.attr("dy", ".35em")
.text(function(d) { return "cell_"+d.cellID; });

// -------------------------------------------------------------------
// update display on multiline
var selected_nodes_data = selected_nodes.data();
//if (selected_nodes_data.length != 0){

console.log(selected_nodes_data);

// clear the old plots
d3.select(".volt").selectAll("."+params.mlp_divID).remove();
d3.select(".calc").selectAll("."+params.mlp_divID).remove();
var mlp_params = {};
mlp_params.divID = params.mlp_divID;
mlp_params.margin =  params.mlp_margin;
mlp_params.width =  params.mlp_width;
mlp_params.height = params.mlp_height;

// update volt-plot
mlp_params.classname =  "volt";
plotMultiLine6(selected_nodes_data, mlp_params);

// update calc-plot
mlp_params.classname =  "calc";
plotMultiLine6(selected_nodes_data, mlp_params);
//}
}
}


//  var brush_elem = petridish.append("g").attr("class", "brush").call(brush);
petridish.classed("brush",true).call(brush);
// =====================================================================================================
// =====================================================================================================

// ---------------------------------------------------------------------------
// Spatial distribution of cells

// associate msrmt_data to visual_elements

var clicked_cells = [];
var cell = petridish.selectAll(".cell").data(cellData);

cell.enter()
.append("svg:g")
.attr("class","g_cell")
  .append("circle") // define attributes of the visual elements being added.
  .attr("class","cell")
  .attr("id",function(d) { return "cell_"+d.cellID; })
    .attr("cx", function(d) { return pd_xScale(d.xPos); })
    .attr("cy", function(d) { return pd_yScale(d.yPos);})
    .attr("pointer-events", "all")
    .on("click", function(d){
                    //alert("Hello! Removed Plots!: enter()");

                    // store new data
                    clicked_cells.push(this.__data__);

                    // erase the existing singlelineplots
                    d3.select(".volt").selectAll("."+params.slp_divID).remove();
                    d3.select(".calc").selectAll("."+params.slp_divID).remove();

                    // update with new single-line plots.
                    var slp_params = {};
                    slp_params.divID = params.slp_divID;
                    slp_params.margin =  params.slp_margin;
                    slp_params.width =  params.slp_width;
                    slp_params.height = params.slp_height;

                    // update volt-singleline
                    slp_params.classname =  "volt";
                    plotMultiLine6(clicked_cells, slp_params);
                    // update calc-singleline
                    slp_params.classname =  "calc";
                    plotMultiLine6(clicked_cells, slp_params);

                    })
    .attr('fill', 'red')
    //.style("fill","steelblue")
    .attr("r", 20)
    .transition(10) // let the circles expand in
      .attr("r",10)
      .attr("fill",function(d){return "url(#gradient_"+pd_classname+"_"+d.cellID+")"; });


cell.exit()
  .transition() // let circles shrink-out
    .attr('r',0)
    .remove();


// update the DATA stored for each circle!

petridish.selectAll(".g_cell").datum(function(d){d.cx = pd_xScale(d.xPos); d.cy = pd_yScale(d.yPos); return d; });
petridish.selectAll("circle").datum(function(d){d.cx = pd_xScale(d.xPos); d.cy = pd_yScale(d.yPos); return d; });
    // --------------------------------------------------------------------
    // cell-edges


var edge = petridish.selectAll("path").data(linkData);


// -----------------------------------------------------------------------------
// representing edges as CURVED-LINES!
    edge.enter()
    .append("path")
       .attr("class", "cell_link")
       .style("stroke-width",function(d){
          var edge_width = Math.ceil(max_edge_thickness*d.value);
          console.log(edge_width);
          return edge_width;
       })
       .attr("d",function(d){
         var src = d.src;
         var src_x =d3.select("."+params.pd_class_name).selectAll("#cell_"+src).attr("cx");
         var src_y =d3.select("."+params.pd_class_name).selectAll("#cell_"+src).attr("cy");
         //console.log(src_x, src_y);
        var dest = d.dest;
         var dest_x =d3.select("."+params.pd_class_name).selectAll("#cell_"+dest).attr("cx");
         var dest_y =d3.select("."+params.pd_class_name).selectAll("#cell_"+dest).attr("cy");
        //console.log(dest_x, dest_y);


        var sum_x = ((+src_x) + (+dest_x) )/2;
        var sum_y = ( (+src_y) + (+dest_y) )/2;

        var diff_x = ( (+dest_x) - (+src_x))/2;
        var diff_y = ((+dest_y) - (+src_y))/2;

        var len_diff = Math.sqrt(diff_x*diff_x + diff_y*diff_y);
        var len_sum = Math.sqrt(sum_x*sum_x + sum_y*sum_y);


        var ctrl_x;
        var ctrl_y;
        //console.log(ctrl_x, ctrl_y);
        switch (edge_type){
              case "arc":
                  ctrl_x = (sum_x) + (diff_x/ len_diff); //
                  ctrl_y = (sum_y) + (diff_y/ len_diff);
                  path_dir = "M" + src_x + "," + src_y + "Q" + ctrl_x + "," + ctrl_y + " " + dest_x + "," + dest_y;
                  break;

              case "line":
                  path_dir = "M" + src_x + " " + src_y + " " + dest_x + " " + dest_y;
                  break;

              case "curved":
                  ctrl_x = Math.abs(diff_x) ;
                  ctrl_y = Math.abs(diff_y);

                  path_dir = "M" + src_x + "," + src_y + "S" + ctrl_x + "," + ctrl_y + " " + dest_x + "," + dest_y;
                  break;
          }

          console.log(path_dir);
          return path_dir;
          });


//==============================================================================
//==============================================================================
// Visual Effects: time-series based color-variation of the Cells/Circles in Petridish

// Methodology: Define <svg_defs:radial_gradiaent> for each cell and modify its properties as per the data.
// Similar to associating each cell to the data, we associate each <svg:defs> to the data.
// Thus, indirectly, we associate each <svg:defs> to each cell.

var cell_grad = petridish.selectAll("defs").data(cellData);
cell_grad.enter()
.append("defs") // cellData[i]: is associated with each <defs>
.append("radialGradient")
.attr("id",function(d){return "gradient_"+pd_classname+"_"+d.cellID;}) // accesses data attached to parent <svg:defs>
.attr("cx", "0.5").attr("cy", "0.5").attr("fx", "0.5").attr("fy", "0.5")
.attr("r","0.5");


var radial_gradients = petridish.selectAll("radialGradient");
switch(params.pd_class_name){
case "volt":
  radial_gradients.append("stop").attr("offset", "0%").style("stop-color", "red");
  break;
case "calc":
  radial_gradients.append("stop").attr("offset", "0%").style("stop-color", "green");
  break;
}
radial_gradients.append("stop").attr("offset", "100%").style("stop-color", "white");

//*/

}
