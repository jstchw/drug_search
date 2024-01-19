import CasesTimeSeriesChart from "../ApexChart/CasesTimeSeriesChart"
import { Nav, OverlayTrigger, Popover} from "react-bootstrap"
import WordCloudChart from "../WordCloudChart/WordCloudChart"
import { ExclamationCircle, QuestionCircle } from "react-bootstrap-icons"
import TermChart from "../ApexChart/TermChart";
import { Carousel } from 'react-responsive-carousel';
import 'react-responsive-carousel/lib/styles/carousel.min.css';
import React from "react";
import {Cloud, List} from "@phosphor-icons/react"
import {useUrlParams} from "../../hooks/useUrlParams";
import {URLParams} from "../../types";

const getChartWarning = (params: URLParams, text: string, status: string) => {
    const popover = (
        <Popover>
            <Popover.Body>
                <div>{text}</div>
            </Popover.Body>
        </Popover>
    )

    return (
        <OverlayTrigger trigger={['hover', 'focus']} placement={'bottom'} overlay={popover}>
            <div style={{cursor: 'default'}}>
                {params.searchBy !== 'side_effects' ?
                    <span>Side effects</span> :
                    <span>Substances</span>
                }
                <sup>
                    {status === 'warn' &&
                        <ExclamationCircle className={'mx-1'}/>}
                    {status === 'info' &&
                        <QuestionCircle className={'mx-1'}/>}

                </sup>
            </div>
        </OverlayTrigger>
    )
}

const ChartSection = () => {
    const adverse_effects_text = 'Correlation does not imply causation.'
    const { params } = useUrlParams()

    const [index, setIndex] = React.useState<number>(0);

    return (
        <div className={'mt-4'}>
            <div>
                <h4 className={'text-center'}>Reports over time</h4>
                <CasesTimeSeriesChart/>
            </div>
            <div className={'mt-4'}>
                <h4 className={'text-center'}>
                    {getChartWarning(params, adverse_effects_text, 'warn')}
                </h4>
            </div>
            <Nav variant="tabs" defaultActiveKey={index} className={'mt-3'}>
                <Nav.Item>
                    <Nav.Link className={'d-flex align-items-center'} eventKey="0" onClick={() => setIndex(0)}>
                        <List weight={'light'}/>
                        <div className={'vr mx-2'}/>
                        Term Chart
                    </Nav.Link>
                </Nav.Item>
                <Nav.Item>
                    <Nav.Link className={'d-flex align-items-center'} eventKey="1" onClick={() => setIndex(1)}>
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
                selectedItem={index}
            >
                <div>
                    <TermChart />
                </div>
                <div className={'h-100'}>
                    <WordCloudChart />
                </div>
            </Carousel>
        </div>
    );
}

export default ChartSection;