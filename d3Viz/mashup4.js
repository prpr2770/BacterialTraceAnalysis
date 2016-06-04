
function mashup4(){
//=============================================================================
//=============================================================================
// Call the functions to plot!

// ----------------------------------------------------------------------------
// multilineplot parameters
var mlp_margin = {top: 20, right: 80, bottom: 30, left: 80},
mlp_width = 1200 - mlp_margin.left - mlp_margin.right,
mlp_height = 400 - mlp_margin.top - mlp_margin.bottom;

// multilineplot
var mlp_divID = "#multilineplot";
plotMultiLine(cellData, mlp_divID, mlp_margin, mlp_width, mlp_height);

// ----------------------------------------------------------------------------
// singlelineplot parameters
var slp_margin = {top: 20, right: 80, bottom: 30, left: 80},
slp_width = 1200 - slp_margin.left - slp_margin.right,
slp_height = 400 - slp_margin.top - slp_margin.bottom;

// singlelineplot
var slp_divID = "#lineplot";
var linePlot_data = cellData[0].voltSig;
var linePlot_name = "cell_"+cellData[0].cellID;
var slp_graph_props = plotSingleLine(linePlot_data, linePlot_name, slp_divID, slp_margin, slp_width, slp_height);

var slp_xScale = slp_graph_props[0];
var slp_yScale = slp_graph_props[1];
var slp_xAxis = slp_graph_props[2];
var slp_yAxis = slp_graph_props[3];
var slp_linePlot = slp_graph_props[4];

//[xScale, yScale, xAxis, yAxisLeft, linePlot] = plotSingleLine(linePlot_data, slp_divID, slp_margin, slp_width, slp_height);

// ----------------------------------------------------------------------------
// Time-Lapse Bar: Indicator

var tl_divID = ".tl_div";
var tl_scale = plotTimeLapseBar(tl_divID, slp_margin, slp_width, slp_height );

// =====================================================================================================
// ---------------------------------------------------------------------------
      // Petridish layout
      var margin = {top: 20, right: 50, bottom: 20, left: 50};

      var width = 1800, height = 800;

      var  petriWidth = width - margin.left - margin.right;
      var petriHeight = height - margin.top - margin.bottom;

      var pd_div_Class = ".petridish";

// plotPetriDish(cellData,pd_divID,pd_margin, pd_width, pd_height);
// function plotPetriDish(cellData,pd_div_Class,pd_margin, pd_width, pd_height){}

var petridish = d3.select(pd_div_Class).attr("width", width ).attr("height", height );
/*
// =====================================================================================================
// =====================================================================================================
// Defining Brush:
var nodes = d3.selectAll(".cell");
var brush = petridish.append("g")
    .attr("class", "brush")
      .call(
        d3.svg.brush()
        .x(d3.scale.identity().domain([0, width]))
        .y(d3.scale.identity().domain([0, height]))
        .on("brush", function() {
          var extent = d3.event.target.extent();

          // change NODE to the CIRCLE_CELL identifier.
          nodes.classed("selected", function(d) {  // change the d.x to d.CX or something, according to CELL info.
            return extent[0][0] <= pd_xScale(d.xPos) && pd_xScale(d.xPos) < extent[1][0] && extent[0][1] <= pd_yScale(d.yPos) && pd_yScale(d.yPos) < extent[1][1];
            //return extent[0][0] <= this.cx && this.cx < extent[1][0] && extent[0][1] <= this.cy && this.cy < extent[1][1];
            });//classed

          })// on
        );//call
// =====================================================================================================
// =====================================================================================================
*/

///*
// =====================================================================================================
// =====================================================================================================
// Defining and Using Brush
var brush = d3.svg.brush()
	.x(d3.scale.identity().domain([0, width]))
	.y(d3.scale.identity().domain([0, height]))
	//.on("brushstart", brushstart)
	.on("brush", brushmove)
	.on("brushend", brushend);

var brushCell;

/*
	// Clear the previously-active brush, if any.
	function brushstart() {
		if (brushCell !== this) {
			d3.select(brushCell).call(brush.clear());
			brushCell = this;
		}
	}
*/
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
      var mlp = d3.selectAll(".multiline").remove();
      plotMultiLine(selected_nodes_data, mlp_divID, mlp_margin, mlp_width, mlp_height);
      //}
    }
	}


//  var brush_elem = petridish.append("g").attr("class", "brush").call(brush);
  petridish.classed("brush",true).call(brush);
// =====================================================================================================
// =====================================================================================================
//*/

// ---------------------------------------------------------------------------
// Spatial distribution of cells

      var pd_xScale = d3.scale.linear().range([margin.left, petriWidth]);
      var pd_yScale = d3.scale.linear().range([petriHeight, margin.top]);

      // set domain by considering dataset.
      pd_xScale.domain([d3.min(cellData, function(d){return d.xPos;}), d3.max(cellData, function(d){return d.xPos;})]);
      pd_yScale.domain([d3.min(cellData, function(d){return d.yPos;}), d3.max(cellData, function(d){return d.yPos;})]);

      // associate msrmt_data to visual_elements

      // MODIFY THIS TO INTRODUCE G!
      var cell = petridish.selectAll(".cell").data(cellData);

/*
      cell //.update()
      		.attr("cx", function(d) { return pd_xScale(d.xPos); })
      		.attr("cy", function(d) { return pd_yScale(d.yPos); })
      		.attr("fill",function(d){return "url(#gradient_"+d.cellID+")"; })
          .on("click", function(d){
                        alert("Hello! Remove Plots!: update()");
                        d3.selectAll(".view_plot").remove();

                        // referencing new data to be plotted
                        var line_plot = d3.select(".aLinePlot");

                        // critical line!! -- where the scaleUpdate occurs.
                        line_plot.data(function(){
                          slp_xScale.domain([0,d.voltSig.length]);
                          slp_yScale.domain([d3.min(d.voltSig), d3.max(d.voltSig)]);
                          // update Axis functions
                          slp_xAxis.scale(slp_xScale);
                          slp_yAxis.scale(slp_yScale);

                          // update linePlot function
                          slp_linePlot.x(function(datum,indx) {  return slp_xScale(indx); }).y(function(datum) { return slp_yScale(datum);})

                          return [d.voltSig];
                        });

                        // call the Axes functions
                        d3.select(".x axis").call(slp_xAxis);
                        d3.select(".y axis").call(slp_yAxis);

                        // call the linePlot function
                        line_plot.append("svg:path")
                          .attr("class","view_plot")
                          .attr("d", function(d){return slp_linePlot(d);});
                        });
*/

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

                          /*
                          // removing old traces
                          d3.selectAll(".view_plot_name").remove();
                          d3.selectAll(".view_plot").remove();
                          */


                        // referencing new data to be plotted
                          var line_plot = d3.select(".aLinePlot");

                          // critical line!! -- where the scaleUpdate occurs.
                          line_plot.data(function(){
                              // update scale
                              slp_xScale.domain([0,d.voltSig.length]);
                              slp_yScale.domain([d3.min(d.voltSig), d3.max(d.voltSig)]);
                              // update axis
                              slp_xAxis.scale(slp_xScale);
                              slp_yAxis.scale(slp_yScale);
                              // update linePlot
                              slp_linePlot.x(function(datum,indx) {return slp_xScale(indx); }).y(function(datum) { return slp_yScale(datum);});

                              return [d.voltSig];
                              });

                              // call the Axes functions
                              d3.select(".x axis").call(slp_xAxis);
                              d3.select(".y axis").call(slp_yAxis);

                              // call linePlot function
                              line_plot.append("svg:path")
                                .attr("class","view_plot")
                                .attr("d", function(d){return slp_linePlot(d);});

                                line_plot.append("text")// title of cell, place it near last datapoint!
                                        .attr("class", "view_plot_name")
                                          .attr("transform", function(d) { var indx_last_element = d.length -1; return "translate(" + slp_xScale(indx_last_element) + "," + slp_yScale(d[indx_last_element]) + ")"; })
                                          .attr("x", 3)
                                          .attr("dy", ".35em")
                                          .text(function() { return "cell_"+d.cellID; });

                          })
      		.attr('fill', 'red')
      		//.style("fill","steelblue")
      		.attr("r", 20)
      		.transition(10) // let the circles expand in
      			.attr("r",10)
      			.attr("fill",function(d){return "url(#gradient_"+d.cellID+")"; });


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
/*
// -----------------------------------------------------------------------------
// representing edges as LINES!
          edge.enter()
          .append("line")
             .attr("class", "cell_link")
             .attr("x1", function(d) {
               var cellID = d.src;
               var xpos = d3.select("#cell_"+cellID).attr("cx");
               var ypos = d3.select("#cell_"+cellID).attr("cy");

               // update y1 for line
               d3.select(this).attr("y1", ypos);
               return xpos;
             })
             .attr("x2", function(d) {
               var cellID = d.dest;
               //obtain position values from attributes of circle, and not the data!
               var xpos = d3.select("#cell_"+cellID).attr("cx");
               var ypos = d3.select("#cell_"+cellID).attr("cy");

               // update y1 for line
               d3.select(this).attr("y2", ypos);
               return xpos;
             })
             .attr("stroke-width",function(d){return d.value; });
*/

///*
// -----------------------------------------------------------------------------
// representing edges as CURVED-LINES!
          edge.enter()
          .append("path")
             .attr("class", "cell_link")
             .attr("d",function(d){
               var src = d.src;
               var src_x =d3.select("#cell_"+src).attr("cx");
               var src_y =d3.select("#cell_"+src).attr("cy");
               console.log(src_x, src_y);
              var dest = d.dest;
               var dest_x =d3.select("#cell_"+dest).attr("cx");
               var dest_y =d3.select("#cell_"+dest).attr("cy");
              console.log(dest_x, dest_y);
               var ctrl_x = Math.abs(src_x - dest_x) /2 ;
               var ctrl_y = Math.abs(src_y - dest_y) /2 ;
               console.log(ctrl_x, ctrl_y);
               return "M" + src_x + "," + src_y + "S" + ctrl_x + "," + ctrl_y + " " + dest_x + "," + dest_y;


             })
             .attr("stroke-width",function(d){return d.value; });
// ---------------------------------------------------------------------------
//*/


//==============================================================================
//==============================================================================
//d3.selectAll(".cell").attr("pointer-events", "all");
      // ---------------------------------------------------------------------------
      // Visual Effects: time-series based color-variation of the Cells/Circles in Petridish

      // Methodology: Define <svg_defs:radial_gradiaent> for each cell and modify its properties as per the data.
      // Similar to associating each cell to the data, we associate each <svg:defs> to the data.
      // Thus, indirectly, we associate each <svg:defs> to each cell.

      var cell_grad = petridish.selectAll("defs").data(cellData);
      cell_grad.enter()
      	.append("defs") // cellData[i]: is associated with each <defs>
      		.append("radialGradient")
      			.attr("id",function(d){return "gradient_"+d.cellID;}) // accesses data attached to parent <svg:defs>
      			.attr("cx", "0.5").attr("cy", "0.5").attr("fx", "0.5").attr("fy", "0.5")
      			.attr("r","0.5");


      var radial_gradients = d3.selectAll("radialGradient");
      radial_gradients.append("stop").attr("offset", "0%").style("stop-color", "red");
      radial_gradients.append("stop").attr("offset", "100%").style("stop-color", "white");




// =====================================================================================================
// Incorporating variation in graphics over time!

var time_index = 0; // keep track of time-index
var tick_duration = 100;

d3.select("#tick_panel").append("p")
.attr("class","tick_panel_value")
.text(tick_duration);

// when the input range changes update the circle
d3.select("#nTicks").on("input", function() {

        // remove last represented value
          d3.selectAll(".tick_panel_value").remove();

        // update variable value
          tick_duration = +this.value;

        // display/print new value
        d3.select("#tick_panel").append("p")
        .attr("class","tick_panel_value")
        .text(tick_duration);

   //interval = window.setInterval(repeatFunctionAtIntervals, window_interval);
   clearInterval(interval);
   interval = window.setInterval(repeatFunctionAtIntervals, tick_duration);
});


// ----------------------------------------------------------------
// ----------------------------------------------------------------
//var window_interval = tick_duration;
var repeatFunctionAtIntervals = function() { // Repeatedly call this function over given time-interval.
          // ===================================================================
          // Visual Effects: time-series based color-variation of cells.
          //tick_duration
          var transition_rate = tick_duration/2;
          //var transition_rate = window_interval/2;
          d3.select(".tl_circle").transition().duration(transition_rate).attr("cx",tl_scale(time_index));

          // ----------------------------------------------------------------------

          var rad_gradients = d3.selectAll("radialGradient");
          rad_gradients.transition().duration(transition_rate) // change the radius of the radial_gradient at each time-step
            .attr("r",
              function(d){
                var pd_rScale = d3.scale.linear().domain([d3.min(d.voltSig), d3.max(d.voltSig)]);
                pd_rScale.range([0, 1]);
                return pd_rScale(d.voltSig[time_index]); }
            );


            // update time-step : such that it runs as an infinte loop
            time_index = (time_index+1) % (cellData[1].voltSig.length);
            // ===================================================================
};

var interval = window.setInterval(repeatFunctionAtIntervals ,tick_duration);
}
