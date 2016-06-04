var testdata = [1, 2, 3, 4, 5];

var new_circles = d3.select(".dots").selectAll("circle").data(testdata);

new_circles.enter().append("circle")
		.style("fill","steelblue")
		.attr("cx", function(d) { return 10*d; })
		.attr("cy", function(d) { return 10*d; })
		.attr("r", 10);


		/*
		data = [1, 2, 3, 4, 5];

		d3.select(".dots").selectAll("circle").data(data)
			.enter()
				.append("circle")
				.style("fill","steelblue")
				.attr("cx", function(d) { return 10*d; })
				.attr("cy", function(d) { return 10*d; })
				.attr("r", 2.5);
		*/

		data = [{x:1,y:1}, {x:2,y:2}, {x:3,y:3}, {x:4,y:4}];

		var dot = d3.select(".dots").selectAll("circle").data(data);

dot.attr("cx", function(d) { return 10*d.x; })
				.attr("cy", function(d) { return 10*d.y; });
dot.enter().append("circle")
				.attr("cx", function(d) { return 10*d.x; })
				.attr("cy", function(d) { return 10*d.y; })
				.attr("r", 0)
				.transition()
					.attr("r",5);


		//  =====================================================================================================
		// test to learn if data can be accessed at multiple levels below it.
		var strataData = [{x:1,y:1},{x:1,y:2},{x:2,y:3},{x:3,y:4}];
		var strata_div = d3.select(".strata").selectAll("div");

		// first level insertion
		var sub_div = strata_div.data(strataData).enter().append("div");
		// second/third level insertion : can access data associated in first-level.
		sub_div.append("p").append("p").text(function(d){return d.x + d.y;});



// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
var brush = d3.svg.brush()
		.x(pd_xScale)
		.y(pd_yScale)
		.on("brushstart", brushstart)
		.on("brush", brushmove)
		.on("brushend", brushend);

		var brushCell;

		// Clear the previously-active brush, if any.
		function brushstart(p) {
			if (brushCell !== this) {
				d3.select(brushCell).call(brush.clear());
				x.domain(domainByTrait[p.x]);
				y.domain(domainByTrait[p.y]);
				brushCell = this;
			}
		}

		// Highlight the selected circles.
		function brushmove(p) {
			var e = brush.extent();
			svg.selectAll("circle").classed("hidden", function(d) {
				return e[0][0] > d[p.x] || d[p.x] > e[1][0] || e[0][1] > d[p.y] || d[p.y] > e[1][1];
			});
		}

		// If the brush is empty, select all circles.
		function brushend() {
			if (brush.empty()) svg.selectAll(".hidden").classed("hidden", false);
		}


//==================================================================

var width = 960,
    height = 500;

var svg = d3.select("body").append("svg")
    .attr("width", width)
    .attr("height", height);

d3.json("graph.json", function(error, graph) {

  graph.links.forEach(function(d) {
    d.source = graph.nodes[d.source];
    d.target = graph.nodes[d.target];
  });

  var link = svg.append("g")
      .attr("class", "link")
    .selectAll("line")
      .data(graph.links)
    .enter().append("line")
      .attr("x1", function(d) { return d.source.x; })
      .attr("y1", function(d) { return d.source.y; })
      .attr("x2", function(d) { return d.target.x; })
      .attr("y2", function(d) { return d.target.y; });

  var node = svg.append("g")
      .attr("class", "node")
    .selectAll("circle")
      .data(graph.nodes)
    .enter().append("circle")
      .attr("r", 4)
      .attr("cx", function(d) { return d.x; })
      .attr("cy", function(d) { return d.y; });

  var brush = svg.append("g")
      .attr("class", "brush")
      	.call(
					d3.svg.brush()
        	.x(d3.scale.identity().domain([0, width]))
        	.y(d3.scale.identity().domain([0, height]))
        	.on("brush", function() {
          	var extent = d3.event.target.extent();
          	node.classed("selected", function(d) {
            	return extent[0][0] <= d.x && d.x < extent[1][0] && extent[0][1] <= d.y && d.y < extent[1][1];
          		});//classed
        		})// on
					);//call
});// json
