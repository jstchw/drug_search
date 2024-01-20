import {useTermData} from "../../hooks/useTermData";
import {ApexOptions} from "apexcharts";
import {ThemeType, URLParams} from "../../types";
import ReactApexChart from "react-apexcharts";
import {ThemeContext} from "../../contexts/ThemeContext";
import React from "react";
import ReactWordcloud, {OptionsProp} from "react-wordcloud";
import {Carousel} from "react-responsive-carousel";
import {Nav, OverlayTrigger, Popover} from "react-bootstrap";
import {Cloud, List} from "@phosphor-icons/react";
import { SealWarning, ChartLine, SmileyNervous } from "@phosphor-icons/react";
import {useUrlParams} from "../../hooks/useUrlParams";

const cloudOptions: OptionsProp  = {
    enableTooltip: true,
    enableOptimizations: true,
    deterministic: true,
    fontFamily: "trebuchet ms",
    fontSizes: [10, 60],
    fontStyle: "normal",
    fontWeight: "normal",
    padding: 1,
    rotations: 0,
    scale: "log",
    spiral: "archimedean",
    transitionDuration: 1000,
}

const TermCarousel = () => {
    const { theme } = React.useContext(ThemeContext)
    const [carouselIndex, setCarouselIndex] = React.useState<number>(0);
    const { params } = useUrlParams()

    const { data, error } = useTermData()

    if (!data || error) {
        return null
    }

    const totalSideEffectCount = data.reduce((acc, obj) => acc + obj.y, 0)

    const chartData = {
        labels: data.map((obj) => obj.x),
        series: [{
            name: 'Events',
            data: data
        }]
    }

    const apexColors = ["#59768A", "#035363", "#32B2BF", "#D5E0BE", "#CE9062", "#E0AB86", "#C7CE8A", "#6EB585", "#325951", "#6F9F9D"]

    const apexChartOptions: ApexOptions = {
        colors: apexColors,
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
            dropShadow: {
                enabled: true,
                top: 1,
                left: 1,
                blur: 2,
                color: theme === 'dark' ? '#000' : '#000',
                opacity: 0.2,
            }
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
        },
    }

    const cloudData = data.map((item) => {
        return {
            text: item.x,
            value: item.y,
        }
    })

    const cloudCallbacks = {
        // Fix this later type
        getWordColor: (word: { text: string; value: number; }) => {
            const maxFrequency = Math.max(...cloudData.map(w => w.value))
            const frequencyRatio = word.value / maxFrequency

            if (theme === 'dark') {
                const darkGray = [172, 181, 189];
                const brightRed = [255, 0, 0];

                const interpolatedColor = darkGray.map((start, i) => {
                    const end = brightRed[i];
                    if (end !== undefined) {
                        return Math.floor(start + frequencyRatio * (end - start));
                    } else {
                        return start;
                    }
                });

                return `rgb(${interpolatedColor.join(', ')})`;
            } else {
                const redComponent = Math.floor(frequencyRatio * 255);
                return `rgb(${redComponent}, 0, 0)`;
            }
        },
        // Fix this later type
        // onWordClick: (word: { text: string; }) => {
        //     if(props.searchOptions.searchBy === searchTypes[2]) {
        //         setSelectedWord(word.text);
        //         setShowDemographicModal(true);
        //     }
        // }
    };

    const getChartWarning = (params: URLParams) => {
        const popover = (
            <Popover>
                <Popover.Body>
                    <div className={'d-flex align-items-center justify-content-center'}>
                        <ChartLine weight={'light'}/>
                        <div className={'vr mx-2'}/>
                        Correlation does not imply causation.
                    </div>
                </Popover.Body>
            </Popover>
        )

        return (
            <OverlayTrigger trigger={['hover', 'focus']} placement={'bottom'} overlay={popover}>
                <div style={{cursor: 'default'}}>
                    {params.searchBy !== 'side_effects' ?
                        <span className={'d-inline-flex align-items-center'}>
                            <SmileyNervous className={'text-secondary'} weight={'light'} />
                            <SealWarning className={'text-secondary'} weight={'light'}/>
                            <div className={'vr mx-2'} />
                            <span>Side Effects</span>
                        </span>
                        :
                        <span>Substances</span>
                    }
                </div>
            </OverlayTrigger>
        )
    }

    return (
        <>
            <h3 className={'d-flex justify-content-center'}>
                {getChartWarning(params)}
            </h3>
            <Nav variant="tabs" defaultActiveKey={carouselIndex} className={'mt-3'}>
                <Nav.Item>
                    <Nav.Link className={'d-flex align-items-center'} eventKey="0" onClick={() => setCarouselIndex(0)}>
                        <List weight={'light'}/>
                        <div className={'vr mx-2'}/>
                        Term Chart
                    </Nav.Link>
                </Nav.Item>
                <Nav.Item>
                    <Nav.Link className={'d-flex align-items-center'} eventKey="1" onClick={() => setCarouselIndex(1)}>
                        <Cloud weight={'light'}/>
                        <div className={'vr mx-2'}/>
                        Word Cloud
                    </Nav.Link>
                </Nav.Item>
            </Nav>
            <Carousel
                showThumbs={false}
                showIndicators={false}
                showArrows={false}
                showStatus={false}
                selectedItem={carouselIndex}
            >
                <ReactApexChart
                    options={apexChartOptions}
                    series={chartData.series}
                    type={apexChartOptions.chart?.type}
                />
                <ReactWordcloud words={cloudData} options={cloudOptions} callbacks={cloudCallbacks} />
            </Carousel>
        </>
    )
}

export default TermCarousel;