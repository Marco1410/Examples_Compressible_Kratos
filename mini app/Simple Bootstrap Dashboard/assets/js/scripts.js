var newUsersChart = new CanvasJS.Chart("newUsersChart", {
  backgroundColor: "transparent",
  animationEnabled: true,
  theme: "light2",
  title: {
    text: "5",
    verticalAlign: "center",
    horizontalAlign: "center",
    dockInsidePlotArea: true
  },
  subtitles: [{
    text: "New Users",
    verticalAlign: "bottom"
  }],
  data: [{
    type: "doughnut",
    startAngle: -90,
    innerRadius: "90%",
    toolTipContent: "New Users: {y}",
    dataPoints: [
      { y: 5, color: "#fc9e24" },
      { y: 45, color: "#efefef", highlightEnabled: false, toolTipContent: null }
    ]
  }]
});
newUsersChart.render();

var salesChart = new CanvasJS.Chart("salesChart", {
  backgroundColor: "transparent",
  animationEnabled: true,
  theme: "light2",
  title: {
    text: "36",
    verticalAlign: "center",
    horizontalAlign: "center",
    dockInsidePlotArea: true
  },
  subtitles: [{
    text: "Sales",
    verticalAlign: "bottom"
  }],
  data: [{
    type: "doughnut",
    startAngle: -90,
    innerRadius: "90%",
    toolTipContent: "Sales: {y}",
    dataPoints: [
      { y: 36, color: "#25b62f" },
      { y: 14, color: "#efefef", highlightEnabled: false, toolTipContent: null }
    ]
  }]
});
salesChart.render();

var subscribersChart = new CanvasJS.Chart("subscribersChart", {
  backgroundColor: "transparent",
  animationEnabled: true,
  theme: "light2",
  title: {
    text: "12",
    verticalAlign: "center",
    horizontalAlign: "center",
    dockInsidePlotArea: true
  },
  subtitles: [{
    text: "Subscribers",
    verticalAlign: "bottom"
  }],
  data: [{
    type: "doughnut",
    startAngle: -90,
    innerRadius: "90%",
    toolTipContent: "Subscribers: {y}",
    dataPoints: [
      { y: 12, color: "#ef5c61" },
      { y: 38, color: "#efefef", highlightEnabled: false, toolTipContent: null }
    ]
  }]
});
subscribersChart.render();

var axisY = {
  lineThickness: 0,
  gridColor: "#def5fc",
  tickLength: 0
};

var incomeChart = new CanvasJS.Chart("incomeChart", {
  backgroundColor: "transparent",
  animationEnabled: true,
  axisX: {
    lineThickness: 0,
    tickThickness: 0
  },
  axisY: axisY,
  data: [{
    color: "#fc9e24",
    dataPoints: [
      { label: "S", y: 50 },
      { label: "M", y: 23 },
      { label: "T", y: 78 },
      { label: "W", y: 35 },
      { label: "T", y: 21 },
      { label: "F", y: 55 },
      { label: "S", y: 25 },
      { label: "S", y: 18 },
      { label: "M", y: 68 },
      { label: "T", y: 82 },
    ]
  }]
});
incomeChart.render();

var userStatisticsChart = new CanvasJS.Chart("userStatisticsChart", {
  backgroundColor: "transparent",
  animationEnabled: true,
  theme: "light2",
  title:{
    text: "New Users",
    fontSize: 18,
    fontColor: "#0dcaf0",
    fontWeight: "normal"
  },
  axisX: {
    lineThickness: 0,
    tickLength: 0,
    labelFormatter: function(e) {
      return "";
    }
  },
  axisY: axisY,
  data: [{
    type: "splineArea", //change type to bar, line, area, pie, etc
    markerSize: 0,
    color: "#0dcaf0",
    fillOpacity: 0.5,
    xValueFormatString: "MMM YYYY",
    dataPoints: [
      { x: new Date(2020, 00, 01), y: 74553 },
      { x: new Date(2020, 01, 01), y: 75374 },
      { x: new Date(2020, 02, 01), y: 76065 },
      { x: new Date(2020, 03, 01), y: 76521 },
      { x: new Date(2020, 04, 01), y: 76695 },
      { x: new Date(2020, 05, 01), y: 76381 },
      { x: new Date(2020, 06, 01), y: 76693 },
      { x: new Date(2020, 07, 01), y: 75624 },
      { x: new Date(2020, 08, 01), y: 76367 },
      { x: new Date(2020, 09, 01), y: 75426 },
      { x: new Date(2020, 10, 01), y: 76024 },
      { x: new Date(2020, 11, 01), y: 75960 }
    ]
  }]
});
userStatisticsChart.render();

var monthlySalesChart = new CanvasJS.Chart("monthlySalesChart", {
  theme: "dark2",
  backgroundColor: "transparent",
  animationEnabled: true,
  title: {
    text: "$890,800",
    horizontalAlign: "left"
  },
  axisX: {
    lineThickness: 0,
    tickLength: 0,
    labelFormatter: function(e) {
      return "";
    }
  },
  axisY: axisY,
  data: [{
    type: "splineArea",
    color: "#ffffff",
    lineColor: "#ffffff",
    fillOpacity: 0.8,
    markerSize: 0,
    xValueFormatString: "MMM YYYY",
    yValueFormatString: "$#,###.##",
    dataPoints: [
      { x: new Date(2020, 00, 01), y: 77000 },
      { x: new Date(2020, 01, 01), y: 76900 },
      { x: new Date(2020, 02, 01), y: 75800 },
      { x: new Date(2020, 03, 01), y: 73000 },
      { x: new Date(2020, 04, 01), y: 70000 },
      { x: new Date(2020, 05, 01), y: 71000 },
      { x: new Date(2020, 06, 01), y: 70500 },
      { x: new Date(2020, 07, 01), y: 73400 },
      { x: new Date(2020, 08, 01), y: 76300 },
      { x: new Date(2020, 09, 01), y: 74700 },
      { x: new Date(2020, 10, 01), y: 75300 },
      { x: new Date(2020, 11, 01), y: 76900 }
    ]
  }]
});
monthlySalesChart.render();


document.getElementById("exportUserChart").addEventListener("click", function() {
  userStatisticsChart.exportChart();
});

document.getElementById("printUserChart").addEventListener("click", function() {
  userStatisticsChart.print();
});