import React from "react";
import ReactApexChart from "react-apexcharts";


// This function sorts terms by popularity and crops the rest
const processDataForChart = (eventDict, totalCount, type) => {
    const maxTerms = 10


    if (type === 'all_groups') {
        let popularTerms = {}
        let sortedDict = Object.entries(eventDict).sort((a, b) => b[1] - a[1])

        for (let [term, count] of sortedDict.slice(0, maxTerms)) {
            popularTerms[term] = count;
        }

        popularTerms = Object.entries(popularTerms).map(([term, count]) => ({
            x: term,
            y: ((count / totalCount) * 100).toPrecision(3),
        }))

        return popularTerms;

    }
    if (type === 'events_over_time') {

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

    if (props.type === 'all_groups') {
        const chartData = {
            labels: processedData.map((x) => x.x),
            series: [{
                name: 'Percentage',
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
                width: '100%',
                toolbar: {
                    show: false,
                }
            },
            plotOptions: {
                bar: {
                    horizontal: false,
                    distributed: true,
                    columnWidth: '95%'
                }
            },
            title : {
                align: 'center',
                text: 'Adverse effects for all groups',
                style: {
                    fontSize: '20px',
                    fontWeight: 'light',
                }
            },
            yaxis: {
                labels: {
                    show: true,
                    formatter: (value) => {
                        return `${value}%`
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
                    type: 'vertical', // Change this to 'horizontal' for horizontal gradient
                    shadeIntensity: 0.4,
                    gradientToColors: undefined,
                    inverseColors: true,
                    opacityFrom: 1,
                    opacityTo: 0.5,
                    stops: [0, 100]
                },
            },
            stroke: {
                show: true,
                width: 3, // This sets the stroke colors
            },
            dataLabels: {
                enabled: false
            },
        }

        return (
            <ReactApexChart
                options={optionsAllGroups}
                series={chartData.series}
                type={optionsAllGroups.chart.type}
            />
        )
    }

    if (props.type === 'events_over_time') {
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
            title : {
                align: 'center',
                text: 'Reports over time',
                style: {
                    fontSize: '20px',
                    fontWeight: 'light',
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