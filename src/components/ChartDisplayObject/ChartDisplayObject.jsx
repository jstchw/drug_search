import ApexChart from "../ApexChart/ApexChart";
import {Alert} from "react-bootstrap";
import React from "react";

const ChartDisplayObject = (props) => {
    return (
        <>
            {props.eventsOverTime && <ApexChart eventDict={props.eventsOverTime} type={'events_over_time'} />}
            <Alert variant={'warning'}>
                <span>Correlation does not imply causation.</span>
            </Alert>
            {props.searchResults && <ApexChart eventDict={props.termCountDict} totalCount={props.totalCount} type={'all_groups'} />}
        </>
    );
}

export default ChartDisplayObject;