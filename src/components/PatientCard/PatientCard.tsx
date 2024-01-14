import React, {useEffect, useRef, useState} from 'react';
import {Toast, ToastContainer, Badge, Col} from 'react-bootstrap';
import {SearchOptions} from "../../types";



const PatientCard = (props: {ADEs: number; searchOptions: SearchOptions}) => {
    const searchOptionsRef = useRef(props.searchOptions);
    const [show, setShow] = useState(false);

    useEffect(() => {
        if(Object.values(props.searchOptions).some(value => value !== undefined))
        setShow(true)
    }, [props.searchOptions]);

    return (
        <ToastContainer
            className={'patient-card m-3'}
            position={'bottom-end'}
            containerPosition={'fixed'}
        >
            <Toast show={show} delay={3000} style={{maxWidth: '200px'}}>
                <Toast.Header className={'justify-content-end'} closeButton={false}>
                    <span>Patient Card</span>
                </Toast.Header>
                <Toast.Body>
                    {/*{Object.keys(searchOptions).map((key, index) => {*/}
                    {/*    return (*/}
                    {/*        <React.Fragment key={index}>*/}
                    {/*            <Col className={'fs-5 d-flex justify-content-start align-items-center'}>*/}
                    {/*                <Badge>{key.charAt(0).toUpperCase() + key.slice(1)}</Badge>*/}
                    {/*                &nbsp;*/}
                    {/*                {searchOptions[key as keyof PatientCardDisplayOptions]}*/}
                    {/*            </Col>*/}
                    {/*        </React.Fragment>*/}
                    {/*    )*/}
                    {/*})}*/}
                    {Object.keys(searchOptionsRef.current).map((key, value) => {
                        return (
                            <React.Fragment key={key}>
                                {value &&
                                    <Col className={'fs-5 d-flex justify-content-start align-items-center'}>
                                        <Badge>{key.charAt(0).toUpperCase() + key.slice(1)}</Badge>
                                        &nbsp;
                                        {value.toString()}
                                    </Col>
                                }
                            </React.Fragment>
                        )}
                    )}
                </Toast.Body>
            </Toast>
        </ToastContainer>
    )
}

export default PatientCard