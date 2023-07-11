import React from "react";
import ReactApexChart from "react-apexcharts";

const chartTypes = [
    'searched_group',
    'events_over_time',
]


// This function sorts terms by popularity and crops the rest
const processDataForChart = (eventDict, totalCount, type) => {
    const maxTerms = 10


    if (type === chartTypes[0]) {
        let popularTerms = {}
        let sortedDict = Object.entries(eventDict).sort((a, b) => b[1] - a[1])

        for (let [term, count] of sortedDict.slice(0, maxTerms)) {
            popularTerms[term] = count;
        }

        popularTerms = Object.entries(popularTerms).map(([term, count]) => ({
            x: term,
            y: count,
        }))

        return popularTerms;

    }
    if (type === chartTypes[1]) {

        let yearlyData = {}
        eventDict.forEach(entry => {
            let year = entry.time.substring(0, 4)
            if (yearlyData[year]) {
                yearlyData[year] += entry.count
            } else {
                yearlyData[year] = entry.count;
            }
        })

        return Object.entries(yearlyData).map(([time, count]) => ({
            x: time,
            y: count
        }))
    }
}


const ApexChart = (props) => {

    const processedData = processDataForChart(props.eventDict, props.totalCount, props.type)

    if (props.type === chartTypes[0]) {
        const chartData = {
            labels: processedData.map((x) => x.x),
            series: [{
                name: 'Events',
                data: processedData
            }]
        }

        const optionsSearchedGroup = {
            labels: chartData.labels,
            legend: {
                show: false,
            },
            chart: {
                type: 'bar',
                width: '100%',
                toolbar: {
                    show: false,
                }
            },
            plotOptions: {
                bar: {
                    horizontal: true,
                    distributed: true,
                    barHeight: '95%',
                }
            },
            yaxis: {
                labels: {
                    show: true,
                    style: {
                        fontSize: '14px',
                    },
                    formatter: (value) => {
                        value = String(value)
                        return value.charAt(0).toUpperCase() + value.slice(1).toLowerCase()
                    }
                }
            },
            xaxis: {
                labels: {
                    show: false
                }
            },
            fill: {
                type: 'gradient',
                gradient: {
                    shade: 'light',
                    type: 'horizontal', // Change this to 'horizontal' for horizontal gradient
                    shadeIntensity: 0.25,
                    gradientToColors: undefined,
                    inverseColors: true,
                    opacityFrom: 1,
                    opacityTo: 0.75,
                    stops: [0, 100]
                },
            },
            stroke: {
                show: true,
                width: 3, // This sets the stroke colors
            },
            dataLabels: {
                enabled: true,
                // Converts the value to a percentage
                formatter: (value) => {
                    return `${((value / props.totalCount) * 100).toPrecision(3)}%`
                }
            },
            tooltip: {
                y: {
                    formatter: (value) => {
                        return `${value.toLocaleString()} from ${props.totalCount.toLocaleString()}`
                    }
                }
            }
        }

        return (
            <ReactApexChart
                options={optionsSearchedGroup}
                series={chartData.series}
                type={optionsSearchedGroup.chart.type}
            />
        )
    } else if (props.type === chartTypes[1]) {
        const chartData = {
            series: [{
                name: 'Reports',
                data: processedData
            }]
        }

        const optionsReportsOverTime = {
            legend: {
                show: true,
                position: 'bottom',
            },
            chart: {
                type: 'area',
                toolbar: {
                    show: false,
                    tools: {
                        zoom: false,
                        zoomin: false,
                        zoomout: false
                    }
                }
            },
            stroke: {
                curve: 'smooth'
            },
            fill: {
                type: 'gradient'
            },
            dataLabels: {
                enabled: false
            },
        }

        return (
            <ReactApexChart
                options={optionsReportsOverTime}
                series={chartData.series}
                type={optionsReportsOverTime.chart.type}
            />
        )
    }
}

export default ApexChart;