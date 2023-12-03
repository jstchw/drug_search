import React from "react";
import ReactApexChart from "react-apexcharts";
import {ThemeContext} from "../../contexts/ThemeContext";
import {ChartDataPoint, Results, ThemeType, TimeEventData} from "../../types";
import {ApexOptions} from "apexcharts";

const chartTypes = [
    'searched_group',
    'events_over_time',
]

const processTermData = (data: Results): ChartDataPoint[] => {
    const entries = Object.entries(data)
    const sortedEntries = entries.sort((a, b) => b[1] - a[1])

    return sortedEntries.slice(0, 10).map(([term, count]) => ({
        x: term,
        y: count,
    }))
}

const processYearData = (data: TimeEventData[]): ChartDataPoint[] => {
    const yearlyData = data.reduce((acc: { [year: string]: number }, entry: TimeEventData) => {
        const year = entry.time.substring(0, 4) // Accessing the year from the time string (e.g. 2019-01-01)
        acc[year] = (acc[year] || 0) + entry.count // If the year exists in the accumulator, add the count to it, otherwise set it to 0 and add the count to it (this is to avoid undefined errors)
        return acc
    }, {})

    const sortedEntries: [string, number][] = Object.entries(yearlyData).sort((a, b) => a[0].localeCompare(b[0]))

    return sortedEntries.map(([time, count]): ChartDataPoint => ({
        x: time,
        y: count,
    }))
}

const ApexChart = (props: { eventDict: TimeEventData[] | Results; totalCount: number; type: string}) => {
    const { theme } = React.useContext(ThemeContext)
    // const processedData = processDataForChart(props.eventDict, props.type)

    if (props.type === chartTypes[0]) {
        const processedData = processTermData(props.eventDict as Results)

        const chartData = {
            labels: processedData.map((x) => x.x),
            series: [{
                name: 'Events',
                data: processedData
            }]
        }

        const options: ApexOptions = {
            colors: ["#59768A", "#035363", "#32B2BF", "#D5E0BE", "#CE9062", "#E0AB86", "#C7CE8A", "#6EB585", "#325951", "#6F9F9D"],
            theme: {
                mode: theme as ThemeType,
            },
            labels: chartData.labels,
            legend: {
                show: false,
            },
            chart: {
                type: "bar",
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
                    formatter: (value: number) => {
                        const stringValue = String(value)
                        return stringValue.charAt(0).toUpperCase() + stringValue.slice(1).toLowerCase()
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
                formatter: (value:number) => {
                    return `${((value / props.totalCount) * 100).toPrecision(3)}%`
                },
            },
            tooltip: {
                y: {
                    formatter: (value: number) => {
                        return `${value.toLocaleString()} from ${props.totalCount.toLocaleString()}`
                    }
                }
            }
        }

        return (
            <ReactApexChart
                options={options}
                series={chartData.series}
                type={options.chart?.type}
            />
        )
    } else if (props.type === chartTypes[1]) {
        const processedData = processYearData(props.eventDict as TimeEventData[])

        const chartData = {
            series: [{
                name: 'Reports',
                data: processedData
            }]
        }

        const options: ApexOptions = {
            theme: {
                mode: theme as ThemeType,
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
                type: 'category',
                labels: {
                    style: {
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
                    formatter: (value) => {
                        return value.toLocaleString('en-US')
                    }
                },
            },
        }

        return (
            // Type specified redundantly because it doesn't work otherwise
            <ReactApexChart
                options={options}
                series={chartData.series}
                type={options.chart?.type}
            />
        )
    } else {
        return <div></div>
    }
}

export default ApexChart;