import ApexChart from "../ApexChart/ApexChart"
import {OverlayTrigger, Popover} from "react-bootstrap"
import WordCloudChart from "../WordCloudChart/WordCloudChart"
import { ExclamationCircle, QuestionCircle } from "react-bootstrap-icons"
import {SearchOptions, Results, TimeEventData} from "../../types";

const getChartWarning = (linkText: string, text: string, status: string) => {
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

const ChartDisplayObject = (props: { eventsOverTime: TimeEventData[]; searchResults: Results; totalCount: number; searchOptions: SearchOptions }) => {
    const adverse_effects_text = 'Correlation does not imply causation.'

    return (
        <div className={'mt-4'}>
            <div>
                <h4>Reports over time</h4>
                {props.eventsOverTime && <ApexChart eventDict={props.eventsOverTime} totalCount={props.totalCount} type={'events_over_time'} />}
            </div>
            <div className={'mt-4'}>
                <h4>{getChartWarning('Adverse effects', adverse_effects_text, 'warn')}</h4>
                {props.searchResults && <ApexChart eventDict={props.searchResults} totalCount={props.totalCount} type={'searched_group'} />}
            </div>
            <div className={'mt-4'}>
                <h4>Word Cloud</h4>
                <WordCloudChart words={props.searchResults} searchOptions={props.searchOptions}/>
            </div>
        </div>
    );
}

export default ChartDisplayObject;