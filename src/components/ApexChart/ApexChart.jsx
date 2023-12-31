import React from "react";
import ReactApexChart from "react-apexcharts";
import {ThemeContext} from "../../contexts/ThemeContext";

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
    const { theme } = React.useContext(ThemeContext)
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
            colors: ["#59768A", "#035363", "#32B2BF", "#D5E0BE", "#CE9062", "#E0AB86", "#C7CE8A", "#6EB585", "#325951", "#6F9F9D"]
            ,
            theme: {
                mode: theme,
            },
            labels: chartData.labels,
            legend: {
                show: false,
            },
            chart: {
                type: 'bar',
                width: '100%',
                toolbar: {
                    show: false,
                },
                background: theme === 'dark' ? '#212529' : '',
            },
            plotOptions: {
                bar: {
                    horizontal: true,
                    distributed: true,
                    barHeight: '100%',
                }
            },
            grid: {
                show: false,
            },
            yaxis: {
                labels: {
                    show: true,
                    style: {
                        fontSize: '14px',
                        colors: theme === 'dark' ? '#ACB5BD' : '',
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
            dataLabels: {
                enabled: true,
                // Converts the value to a percentage
                formatter: (value) => {
                    return `${((value / props.totalCount) * 100).toPrecision(3)}%`
                },
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
            theme: {
                mode: theme
            },
            colors: ['#59768A'],
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
                },
                background: theme === 'dark' ? '#212529' : '',
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
            xaxis: {
                type: 'datetime',
                labels: {
                    style: {
                        theme: 'dark',
                        fontSize: '14px',
                        colors: theme === 'dark' ? '#ACB5BD' : '',
                    },
                },
                tooltip: {
                    enabled: false,
                }
            },
            yaxis: {
                labels: {
                    style: {
                        fontSize: '14px',
                        colors: theme === 'dark' ? '#ACB5BD' : '',
                    },
                },
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