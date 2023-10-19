import ApexChart from "../ApexChart/ApexChart"
import {OverlayTrigger, Popover} from "react-bootstrap"
import React from "react"
import WordCloudChart from "../WordCloudChart/WordCloudChart"
import { ExclamationCircle, QuestionCircle } from "react-bootstrap-icons"
import CausalInterface from "../CausalInterface/CausalInterface";

const getChartWarning = (linkText, text, status) => {
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
                <span>{linkText}</span>
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

const ChartDisplayObject = (props) => {

    const causal_inf_text = 'Causal inference on this website provides insights into potential links between certain' +
        ' drugs and adverse effects. The Z-score alongside each drug indicates the confidence level of any potential' +
        ' connection. A higher Z-score suggests a stronger link, but always consult a healthcare professional for' +
        ' definitive information and guidance.'

    const adverse_effects_text = 'Correlation does not imply causation.'

    return (
        <div className={'mt-4'}>
            <div>
                <h4>Reports over time</h4>
                {props.eventsOverTime && <ApexChart eventDict={props.eventsOverTime} type={'events_over_time'} />}
            </div>
            <div className={'mt-4'}>
                <h4>{getChartWarning('Causal Inference', causal_inf_text, 'info')}</h4>
                <CausalInterface></CausalInterface>
            </div>
            <div className={'mt-4'}>
                <h4>{getChartWarning('Adverse effects', adverse_effects_text, 'warn')}</h4>
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