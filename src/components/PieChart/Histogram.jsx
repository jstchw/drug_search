import React from "react";
import Chart from "react-apexcharts";

// This function sorts terms by popularity and crops the rest
const processDataForChart = (termCountDict) => {
    const sortedTerms = Object.entries(termCountDict).sort((a, b) => b[1] - a[1])

    let popularTerms = {}
    const maxTerms = 10
    let topCount = 0

    for (let [term, count] of sortedTerms.slice(0, maxTerms)) {
        popularTerms[term] = count;
        topCount += count;
    }

    console.log(popularTerms)

    popularTerms = Object.entries(popularTerms).map(([term, count]) => ({
        x: term,
        y: count,
        //y: (count / topCount) * 100
    }))

    //console.log(popularTerms)
    return popularTerms;
}


const Histogram = ({ termCountDict, totalCount }) => {
    const processedData = processDataForChart(termCountDict)

    const chartData = {
        labels: processedData.map((x) => x.x),
        series: [{
            data: processedData
        }]
    }

    const options = {
        labels: chartData.labels,
        legend: {
            show: true,
            position: 'bottom',
        },
        chart: {
            type: 'bar',
        },
        plotOptions: {
            bar: {
                horizontal: true
            }
        },
        title : {
            text: 'Adverse effects for all groups',
        },
        theme: {
        }
    }

    return (
        <Chart
            options={options}
            series={chartData.series}
            type={options.chart.type}
        />
    )
}

export default Histogram;