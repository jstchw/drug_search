import React from "react";
import Chart from "react-apexcharts";

const processDataForChart = (termCountDict, totalCount) => {
    const thresh = 0.7 * totalCount
    const sortedTerms = Object.entries(termCountDict).sort((a, b) => b[1] - a[1])

    let popularTerms = {}
    let otherTermsCount = 0
    let currentCount = 0

    for (const [term, count] of sortedTerms) {
        currentCount += count
        if (currentCount <= thresh) {
            popularTerms[term] = count
        } else {
            otherTermsCount += count
        }
    }

    popularTerms['OTHER'] = otherTermsCount

    return popularTerms
}

const PieChart = ({ termCountDict, totalCount }) => {
    const processedData = processDataForChart(termCountDict, totalCount)

    const chartData = {
        labels: Object.keys(processedData),
        series: Object.values(processedData).map((count) => count / totalCount * 100)
    }

    const options = {
        labels: chartData.labels,
        legend: {
            show: false,
            position: 'bottom',
        },
        chart: {
            type: 'donut',
        },
        tooltip: {
            y: {
                formatter: value => `${value.toFixed(2)}%`
            }
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

export default PieChart;