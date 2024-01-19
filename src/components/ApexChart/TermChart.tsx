import {useTermData} from "../../hooks/useTermData";
import {ApexOptions} from "apexcharts";
import {ThemeType} from "../../types";
import ReactApexChart from "react-apexcharts";
import {ThemeContext} from "../../contexts/ThemeContext";
import React from "react";

const TermChart = () => {
    const { theme } = React.useContext(ThemeContext)
    const { data } = useTermData()
    if (!data) {
        return null
    }
    const totalSideEffectCount = data.reduce((acc, obj) => acc + obj.y, 0)

    const chartData = {
        labels: data?.map((obj) => obj.x),
        series: [{
            name: 'Events',
            data: data
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
                return `${((value / totalSideEffectCount) * 100).toPrecision(3)}%`
            },
        },
        tooltip: {
            y: {
                formatter: (value: number) => {
                    return `${value.toLocaleString()} from ${totalSideEffectCount.toLocaleString()}`
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
}

export default TermChart;