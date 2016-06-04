
function mashup6(cellData,voltLinkData, caLinkData){
//=============================================================================
// introduce multiple-colored line-plots in the singlelineplot div!
//=============================================================================
// Call the functions to plot!

// ----------------------------------------------------------------------------
// multilineplot parameters
var slp_margin = {top: 20, right: 80, bottom: 30, left: 80},
slp_width = 1200 - slp_margin.left - slp_margin.right,
slp_height = 400 - slp_margin.top - slp_margin.bottom;


// multilineplot parameters
var mlp_margin = {top: 20, right: 80, bottom: 30, left: 80},
mlp_width = 1200 - mlp_margin.left - mlp_margin.right,
mlp_height = 400 - mlp_margin.top - mlp_margin.bottom;

// singlelineplot/multilineplot parameters
var mlp_divID = "multiline";
var slp_divID = "singleline";
var linePlot_data = [cellData[0]];

var slp_params = {};
slp_params.divID = slp_divID;
slp_params.margin = slp_margin;
slp_params.width = slp_width;
slp_params.height = slp_height;

var mlp_params = {};
mlp_params.divID = mlp_divID;
mlp_params.margin = mlp_margin;
mlp_params.width = mlp_width;
mlp_params.height = mlp_height;
// ----------------------------------------------------------------------------
// voltSignals
var className = "volt";
// multiline
mlp_params.classname = className;
plotMultiLine6(cellData.slice(1,6), mlp_params);
//singleline
slp_params.classname = className;
plotMultiLine6(linePlot_data, slp_params);
// ----------------------------------------------------------------------------


// caSignals
className = "calc";
// multiline
mlp_params.classname = className;
plotMultiLine6(cellData.slice(1,6), mlp_params);
//singleline
slp_params.classname = className;
plotMultiLine6(linePlot_data, slp_params);


// ----------------------------------------------------------------------------

// Time-Lapse Bar: Indicator

var tl_divID = "tl_div";
className = "volt";
var volt_tl_scale = plotTimeLapseBar6(tl_divID, className, mlp_margin, mlp_width, mlp_height );

className = "calc";
var ca_tl_scale = plotTimeLapseBar6(tl_divID, className, mlp_margin, mlp_width, mlp_height );



// =====================================================================================================
// ---------------------------------------------------------------------------
      // Petridish layout
      var margin = {top: 20, right: 50, bottom: 20, left: 50};

      var width = 1400, height = 800;

      var pd_divID = "petridish";

//var divIDs = [slp_divID, mlp_divID, pd_divID];

// parameters to be sent into plotPetridish

var params = {};
params.pd_width = width;
params.pd_height = height;
params.pd_margin = margin;
params.pd_divID = pd_divID;

params.slp_divID = slp_divID;
params.slp_width = slp_width;
params.slp_height = slp_height;
params.slp_margin = slp_margin;

params.mlp_divID = mlp_divID;
params.mlp_width = mlp_width;
params.mlp_height = mlp_height;
params.mlp_margin = mlp_margin;


// voltage petridish
var pd_class_name = "volt";
params.pd_class_name = pd_class_name;
plotPetridish(voltLinkData, cellData, params);



// calcium petridish
pd_class_name = "calc";
params.pd_class_name = pd_class_name;
plotPetridish(caLinkData, cellData, params);


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
var max_time_index = Math.min(cellData[1].voltSig.length, cellData[1].caSig.length);
var repeatFunctionAtIntervals = function() { // Repeatedly call this function over given time-interval.
      // ===================================================================
      // Visual Effects: time-series based color-variation of cells.
      //tick_duration
      var transition_rate = tick_duration/2;
      //var transition_rate = window_interval/2;
      d3.select(".volt").selectAll(".tl_circle").transition().duration(transition_rate).attr("cx",volt_tl_scale(time_index));
      d3.select(".calc").selectAll(".tl_circle").transition().duration(transition_rate).attr("cx",ca_tl_scale(time_index));
      // ----------------------------------------------------------------------

      var volt_rad_gradients = d3.select(".volt").selectAll("radialGradient");
      volt_rad_gradients.transition().duration(transition_rate) // change the radius of the radial_gradient at each time-step
        .attr("r",
          function(d){
            var pd_rScale = d3.scale.linear().domain([d3.min(d.voltSig), d3.max(d.voltSig)]);
            pd_rScale.range([0, 1]);
            return pd_rScale(d.voltSig[time_index]);
            });


            var ca_rad_gradients = d3.select(".calc").selectAll("radialGradient");
            ca_rad_gradients.transition().duration(transition_rate) // change the radius of the radial_gradient at each time-step
              .attr("r",
                function(d){
                  var pd_rScale = d3.scale.linear().domain([d3.min(d.caSig), d3.max(d.caSig)]);
                  pd_rScale.range([0, 1]);
                  return pd_rScale(d.caSig[time_index]);
                  });

        // update time-step : such that it runs as an infinte loop
        time_index = (time_index+1) % (max_time_index);
        // ===================================================================
      };

var interval = window.setInterval(repeatFunctionAtIntervals ,tick_duration);



}
