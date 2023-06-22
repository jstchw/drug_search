import ApexChart from "../ApexChart/ApexChart";
import {Alert, Col, Container} from "react-bootstrap";
import React from "react";
import WordCloudChart from "../WordCloudChart/WordCloudChart";

const ChartDisplayObject = (props) => {
    return (
        <div className={'mt-4'}>
            {props.eventsOverTime && <ApexChart eventDict={props.eventsOverTime} type={'events_over_time'} />}
            <Alert variant={'warning'}>
                <span>Correlation does not imply causation.</span>
            </Alert>
                {props.searchResults && <ApexChart eventDict={props.termCountDict} totalCount={props.totalCount} type={'searched_group'} />}
                <div>
                    <WordCloudChart ADEArray={props.termCountDict}/>
                </div>
        </div>
    );
}

export default ChartDisplayObject;