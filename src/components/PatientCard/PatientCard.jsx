import React, {useEffect, useRef, useState} from 'react';
import {Toast, ToastContainer, Badge, Row, Col} from 'react-bootstrap';
import {searchSex, searchTypes} from "../OptionModal/OptionModal";


const processSearchOptions = (searchOptions, setSearchOptions, totalADE, searchBy) => {

    let newSearchOptions = {...searchOptions}

    if(searchOptions.age && Array.isArray(searchOptions.age) && searchOptions.age.length === 2 && searchOptions.age[0] && searchOptions.age[1]) {
        newSearchOptions = {
            ...newSearchOptions,
            age: searchOptions.age.join(' - ')
        }
    }
    // If value from searchSex is equal to the value of searchOptions.sex then return the label
    if(searchOptions.sex) {
        const foundSex = searchSex.find(sex => sex.value === searchOptions.sex);
        if (foundSex) {
            newSearchOptions = {
                ...newSearchOptions,
                sex: foundSex.label
            }
        }
    }

    if(searchTypes[2].value !== searchBy) {
        newSearchOptions = {
            ADEs: parseInt(totalADE).toLocaleString('en'),
            ...newSearchOptions,
        }
    } else if(searchTypes[2].value === searchBy) {
        newSearchOptions = {
            ADEs: undefined,
            ...newSearchOptions,
        }
    }
    setSearchOptions(newSearchOptions)
}

const PatientCard = (props) => {
    const { searchBy, ...rest } = props.searchOptions;
    const searchByRef = useRef(searchBy)
    const [searchOptions, setSearchOptions] = React.useState(rest)

    useEffect(() => {
        processSearchOptions(searchOptions, setSearchOptions, props.totalADE, searchByRef.current)
    }, [props.searchOptions, props.totalADE, searchBy])

    const [show, setShow] = useState(false);

    useEffect(() => {
        if(Object.values(searchOptions).some(value => value !== undefined))
        setShow(true)
    }, [searchOptions]);

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