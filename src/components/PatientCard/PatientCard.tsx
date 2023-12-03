import React, {useEffect, useRef, useState} from 'react';
import {Toast, ToastContainer, Badge, Col} from 'react-bootstrap';
import {searchSex, searchTypes} from "../OptionModal/OptionModal";
import {SearchOptions} from "../../types";


type AdditionalOptions = {
    ADEs: string | undefined;
}
// The new type will be the same as SearchOptions but without the searchBy property
type PatientCardDisplayOptions = Omit<SearchOptions, 'searchBy'> & AdditionalOptions;

const processSearchOptions = (searchOptions: PatientCardDisplayOptions, totalADECount: number, searchBy: string) => {

    let newSearchOptions = {...searchOptions}

    if(searchOptions.age && Array.isArray(searchOptions.age) && searchOptions.age.length === 2 && searchOptions.age[0] && searchOptions.age[1]) {
        newSearchOptions = {
            ...newSearchOptions,
            age: {
                value: searchOptions.age.join(' - '),
                index: 0,
                type: ''
            }
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

    if(searchBy !== searchTypes[2]?.value) {
        newSearchOptions = {
            ...newSearchOptions,
            ADEs: totalADECount.toLocaleString('en-US')
        }
    } else if(searchBy === searchTypes[2].value) {
        newSearchOptions = {
            ...newSearchOptions,
            ADEs: undefined
        }
    }

    return newSearchOptions;
}

const PatientCard = (props: {ADEs: number; searchOptions: SearchOptions}) => {
    const { searchBy, ...rest } = props.searchOptions;

    const initialSearchOptionsState: PatientCardDisplayOptions = {
        ...rest,
        ADEs: undefined
    }

    const searchByRef = useRef<string>(searchBy)
    const [searchOptions, setSearchOptions] = React.useState<PatientCardDisplayOptions>(initialSearchOptionsState)

    useEffect(() => {
        const newSearchOptions = processSearchOptions(searchOptions, props.ADEs, searchByRef.current)

        if (JSON.stringify(newSearchOptions) !== JSON.stringify(searchOptions)) {
            setSearchOptions(newSearchOptions);
        }
    }, [searchOptions, props.ADEs]);

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

                        /* ABSOLUTELY REDO THIS */
                        const value = searchOptions[key as keyof PatientCardDisplayOptions];
                        let displayValue;
                        if (typeof value === 'object') {
                            displayValue = value.value;
                        } else {
                            displayValue = value;
                        }
                        if (value !== undefined) {
                            return (
                                <React.Fragment key={index}>
                                    <Col className={'fs-5 d-flex justify-content-start align-items-center'}>
                                        <Badge>{key.charAt(0).toUpperCase() + key.slice(1)}</Badge>
                                        &nbsp;
                                        {displayValue}
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