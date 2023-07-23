import ApexChart from "../ApexChart/ApexChart"
import {OverlayTrigger, Popover} from "react-bootstrap"
import React from "react"
import WordCloudChart from "../WordCloudChart/WordCloudChart"
import { ExclamationCircle } from "react-bootstrap-icons"

const getChartWarning = (linkText) => {
    const popover = (
        <Popover>
            <Popover.Body>
                <div>Correlation does not imply causation.</div>
            </Popover.Body>
        </Popover>
    )

    return (
        <OverlayTrigger trigger={['hover', 'focus']} placement={'bottom'} overlay={popover}>
            <div style={{cursor: 'default'}}>
                <span>{linkText}</span>
                <sup><ExclamationCircle className={'mx-1'}/></sup>
            </div>
        </OverlayTrigger>
    )
}

const ChartDisplayObject = (props) => {
    return (
        <div className={'mt-4'}>
            <div>
                <h4>Reports over time</h4>
                {props.eventsOverTime && <ApexChart eventDict={props.eventsOverTime} type={'events_over_time'} />}
            </div>
            <div className={'mt-4'}>
                <h4>{getChartWarning('Adverse effects')}</h4>
                {props.searchResults && <ApexChart eventDict={props.termCountDict} totalCount={props.totalCount} type={'searched_group'} />}
            </div>
            <div className={'mt-4'}>
                <h4>Word Cloud</h4>
                <WordCloudChart ADEArray={props.termCountDict} searchOptions={props.searchOptions}/>
            </div>
        </div>
    );
}

export default ChartDisplayObject;