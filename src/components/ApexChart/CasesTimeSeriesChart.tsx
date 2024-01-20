import React from "react";
import {ThemeContext} from "../../contexts/ThemeContext";
import {ThemeType} from "../../types";
// import {ApexOptions} from "apexcharts";
// import ReactApexChart from "react-apexcharts";
import {useTimeSeriesData} from "../../hooks/useTimeSeriesData";
import {ApexOptions} from "apexcharts";
import ReactApexChart from "react-apexcharts";
import { Clock } from "@phosphor-icons/react";


const CasesTimeSeriesChart = () => {
    const { theme } = React.useContext(ThemeContext)

    const { data, error } = useTimeSeriesData()

    if (!data || error) {
        return null
    }

    const chartData = {
        series: [{
            name: 'Reports',
            data: data
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
        <>
            <h3 className={'d-flex justify-content-center align-items-center'}>
                <Clock weight={'light'} />
                <div className={'vr mx-2'} />
                Reports over time
            </h3>
            <ReactApexChart
                options={options}
                series={chartData.series}
                type={options.chart?.type}
            />
        </>
    )

}

export default CasesTimeSeriesChart;