import React from "react";
import Chart from "react-apexcharts";

// This function sorts terms by popularity and crops the rest
const processDataForChart = (termCountDict) => {
    const sortedTerms = Object.entries(termCountDict).sort((a, b) => b[1] - a[1])

    let popularTerms = {}
    const maxTerms = 10
    const totalSum = sortedTerms.reduce((sum, term) => sum + term[1], 0);

    for (let [term, count] of sortedTerms.slice(0, maxTerms)) {
        popularTerms[term] = count;
    }



    popularTerms = Object.entries(popularTerms).map(([term, count]) => ({
        x: term,
        y: ((count / totalSum) * 100).toPrecision(3),
    }))

    //console.log(popularTerms)
    return popularTerms;
}


const Histogram = ({ termCountDict, totalCount, type }) => {
    const processedData = processDataForChart(termCountDict)

    const chartData = {
        labels: processedData.map((x) => x.x),
        series: [{
            data: processedData
        }]
    }

    const optionsAllGroups = {
        labels: chartData.labels,
        legend: {
            show: true,
            position: 'bottom',
        },
        chart: {
            type: 'bar',
            width: '100%'
        },
        plotOptions: {
            bar: {
                horizontal: true
            }
        },
        title : {
            text: 'Adverse effects for all groups',
        },
        yaxis: {
            labels: {
                formatter: function (value) {
                    return value + "%";
                }
            }
        }
    }

    const optionsReportsOverTime = {

    }

    return (
        <Chart
            options={optionsAllGroups}
            series={chartData.series}
            type={optionsAllGroups.chart.type}
        />
    )
}

export default Histogram;