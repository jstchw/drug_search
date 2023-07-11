import ApexChart from "../ApexChart/ApexChart";
import {Alert} from "react-bootstrap";
import React from "react";
import WordCloudChart from "../WordCloudChart/WordCloudChart";

const ChartDisplayObject = (props) => {
    return (
        <div className={'mt-4'}>
            <div>
                <h4>Reports over time</h4>
                {props.eventsOverTime && <ApexChart eventDict={props.eventsOverTime} type={'events_over_time'} />}
            </div>
            <Alert variant={'warning'}>
                <span>Correlation does not imply causation.</span>
            </Alert>
            <div>
                <h4>Adverse effects</h4>
                {props.searchResults && <ApexChart eventDict={props.termCountDict} totalCount={props.totalCount} type={'searched_group'} />}
            </div>
            <div>
                <h4>Word Cloud</h4>
                <WordCloudChart ADEArray={props.termCountDict} searchOptions={props.searchOptions}/>
            </div>
        </div>
    );
}

export default ChartDisplayObject;