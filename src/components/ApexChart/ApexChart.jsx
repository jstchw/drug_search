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

        console.log(yearlyData)

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
            },
            plotOptions: {
                bar: {
                    horizontal: false,
                    distributed: true,
                }
            },
            title : {
                text: 'Adverse effects for all groups',
            },
            yaxis: {
                labels: {
                    show: true
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
                    type: 'horizontal',
                    shadeIntensity: 0.25,
                    inverseColors: true,
                    opacityFrom: 0.85,
                    opacityTo: 0.85,
                    stops: [50, 0, 100]
                },
            },
            colors: ['#2E93fA', '#66DA26', '#546E7A', '#E91E63', '#FF9800'],
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
                data: processedData
            }]
        }

        console.log(processedData)

        const optionsReportsOverTime = {
            legend: {
                show: true,
                position: 'bottom',
            },
            chart: {
                type: 'area',
                toolbar: {
                    tools: {
                        zoom: false,
                        zoomin: false,
                        zoomout: false
                    }
                }
            },
            title : {
                text: 'Reports over time',
            },
            stroke: {
                curve: 'smooth'
            },
            fill: {
                type: 'gradient'
            },
            colors: ['#2E93fA', '#66DA26', '#546E7A', '#E91E63', '#FF9800'],
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