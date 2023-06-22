import React, {useEffect, useRef, useState} from 'react';
import {Toast, ToastContainer, Badge, Row, Col} from 'react-bootstrap';
import { searchSex } from "../OptionModal/OptionModal";


const processSearchOptions = (searchOptions, setSearchOptions) => {

    let newSearchOptions = {...searchOptions}

    if(searchOptions.age && searchOptions.age[0] && searchOptions.age[1]) {
        newSearchOptions = {
            ...newSearchOptions,
            age: searchOptions.age.join(' - ')
        }
    }
    // If value from searchSex is equal to the value of searchOptions.sex then return the label
    if(searchOptions.sex) {
        newSearchOptions = {
            ...newSearchOptions,
            sex: searchSex.find(sex => sex.value === searchOptions.sex).label
        }
    }

    setSearchOptions(newSearchOptions)
}

const PatientCard = (props) => {
    const { searchBy, ...rest } = props.searchOptions;
    const [searchOptions, setSearchOptions] = React.useState(rest)

    const [show, setShow] = useState(false);

    useEffect(() => {
        processSearchOptions(searchOptions, setSearchOptions)
        setShow(true)
    }, []);

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
                    <Col className={'fs-5 d-flex justify-content-start align-items-center'}>
                        <Badge>ADEs</Badge>
                        &nbsp;
                        {parseInt(props.totalADE).toLocaleString('en')}
                    </Col>
                    {Object.keys(searchOptions).map((key, index) => {
                        if (searchOptions[key]) {
                            return (
                                <React.Fragment key={index}>
                                    <Col className={'fs-5 d-flex justify-content-start align-items-center'}>
                                        <Badge>{key.charAt(0).toUpperCase() + key.slice(1)}</Badge>
                                        &nbsp;
                                        {searchOptions[key]}
                                    </Col>
                                </React.Fragment>
                            )
                        }
                        return null; // Return null if condition is not met
                    })}
                </Toast.Body>
            </Toast>
        </ToastContainer>
    )
}

export default PatientCard