import React, {useEffect} from 'react';
import {Modal, Form, InputGroup, ToggleButton} from 'react-bootstrap';


// Search types that can be selected in the popover
export const searchTypes = [
    {
        value: 'patient.drug.activesubstance.activesubstancename',
        index: 0,
        label: 'Active substance',
        type: 'searchBy'
    },
    {
        value: 'patient.drug.openfda.brand_name',
        index: 1,
        label: 'Brand Name',
        type: 'searchBy'
    },
    {
        value: 'patient.reaction.reactionmeddrapt.exact',
        index: 2,
        label: 'Side Effect',
        type: 'searchBy'
    },
]

export const searchSex = [
    {
        value: 'patient.patientsex:1',
        index: 0,
        label: 'Male',
        type: 'sex'
    },
    {
        value: 'patient.patientsex:2',
        index: 1,
        label: 'Female',
        type: 'sex'
    }
]

export const searchAgeRange = [
    {
        value: '',
        index: 0,
        type: 'age',
    },
    {
        value: '',
        index: 1,
        type: 'age',
    }
]

export const searchCountry = {
    value: 'occurcountry.exact',
    index: 0,
    label: 'Country',
    type: 'country'
}

const OptionModal = (props) => {
    const [isSexDisabled, setIsSexDisabled] = React.useState(true)
    const [lastSelectedSex, setLastSelectedSex] = React.useState(searchSex[0].value)

    const [isAgeDisabled, setIsAgeDisabled] = React.useState(true)
    const [ageRange, setAgeRange] = React.useState([searchAgeRange[0].value, searchAgeRange[1].value])
    const [lastSelectedAge, setLastSelectedAge] = React.useState([searchAgeRange[0].value, searchAgeRange[1].value])

    // Clear the age range when the age is disabled
    useEffect(() => {
        if(isAgeDisabled) {
            props.setSearchOptions({
                ...props.searchOptions,
                age: undefined
            });
        } else {
            props.setSearchOptions({
                ...props.searchOptions,
                age: lastSelectedAge
            });
        }
    }, [isAgeDisabled]);


    useEffect(() => {
        if(isSexDisabled) {
            props.setSearchOptions({
                ...props.searchOptions,
                sex: undefined
            })
        } else {
            props.setSearchOptions({
                ...props.searchOptions,
                sex: lastSelectedSex
            })
        }
    }, [isSexDisabled])


    const handleOptionsChange = (index, type, value) => {
        switch (type) {
            case 'searchBy':
                props.setSearchOptions({
                    ...props.searchOptions,
                    [type]: searchTypes[index].value
                })

                props.setSelectedSearchOptionIndex({
                    ...props.selectedSearchOptionIndex,
                    [type]: index
                })


                break
            case 'sex':
                setLastSelectedSex(searchSex[index].value)
                props.setSearchOptions({
                    ...props.searchOptions,
                    [type]: searchSex[index].value
                })

                props.setSelectedSearchOptionIndex({
                    ...props.selectedSearchOptionIndex,
                    [type]: index
                })
                break
            case 'age':
                const newAgeRange = [...ageRange]
                newAgeRange[index] = value
                setAgeRange(newAgeRange)
                if(!isAgeDisabled) {
                    setLastSelectedAge(newAgeRange)
                    props.setSearchOptions({
                        ...props.searchOptions,
                        [type]: newAgeRange
                    })
                }
                break
        }
    }

    return (
        <>
            <Modal centered show={props.show} onHide={props.handleClose}>
                <Modal.Header closeButton>
                    <Modal.Title>Options</Modal.Title>
                </Modal.Header>
                <Modal.Body>
                    <Form>
                        <Form.Group className="mb-3 d-flex align-items-center">
                            <div className={'d-flex align-items-center'}>
                                <ToggleButton
                                    type="checkbox"
                                    variant="outline-primary"
                                    checked={true}
                                    value="1"
                                    disabled={true}
                                >
                                    Search by
                                </ToggleButton>
                            </div>
                            <InputGroup className={'flex-grow-1 mx-3'} style={{width: 'auto'}}>
                                <Form.Select
                                    onChange={(e) => handleOptionsChange(e.target.value, 'searchBy')}
                                    value={props.selectedSearchOptionIndex.searchBy}
                                    style={{width: 'auto'}}
                                >
                                    {searchTypes.map((searchType, index) => (
                                        <option
                                            key={index}
                                            value={searchType.index}
                                        >
                                            {searchType.label}
                                        </option>
                                    ))}
                                </Form.Select>
                            </InputGroup>
                        </Form.Group>

                        <Form.Group className="mb-3 d-flex align-items-center">
                            <div className={'d-flex align-items-center'}>
                                <ToggleButton
                                    type="checkbox"
                                    variant="outline-primary"
                                    checked={!isSexDisabled}
                                    value="1"
                                    onClick={(e) => setIsSexDisabled(!isSexDisabled)}
                                >
                                    Sex
                                </ToggleButton>
                            </div>
                            <InputGroup className={'mx-3 flex-grow-1'}>
                                <Form.Select
                                    onChange={e => {handleOptionsChange(e.target.value, 'sex')}}
                                    value={props.selectedSearchOptionIndex.sex}
                                    disabled={isSexDisabled}
                                >
                                    {searchSex.map((sex, index) => (
                                            <option
                                                key={index}
                                                value={sex.index}
                                            >
                                                {sex.label}
                                            </option>
                                ))}
                                </Form.Select>
                            </InputGroup>
                        </Form.Group>

                        <Form.Group className="mb-3 d-flex align-items-center">
                            <div className={'d-flex align-items-center'}>
                                <ToggleButton
                                    type="checkbox"
                                    variant="outline-primary"
                                    checked={!isAgeDisabled}
                                    value="1"
                                    onClick={(e) => setIsAgeDisabled(!isAgeDisabled)}
                                >
                                Age
                                </ToggleButton>
                            </div>

                            <InputGroup className={'mx-3 flex-grow-1'}>
                                <Form.Control
                                    type="number"
                                    placeholder="Min"
                                    value={ageRange[0]}
                                    onChange={(e) => handleOptionsChange(0, 'age', e.target.value)}
                                    disabled={isAgeDisabled}
                                />
                                <InputGroup.Text>-</InputGroup.Text>
                                <Form.Control
                                    type="number"
                                    placeholder="Max"
                                    value={ageRange[1]}
                                    onChange={(e) => handleOptionsChange(1, 'age', e.target.value)}
                                    disabled={isAgeDisabled}
                                />
                            </InputGroup>
                        </Form.Group>
                        <Form.Group>
                            <div className={'d-flex align-items-center'}>
                                <ToggleButton
                                    type="checkbox"
                                    variant="outline-primary"
                                    value="1"
                                >
                                    Country
                                </ToggleButton>
                                <InputGroup className={'mx-3 flex-grow-1'}>
                                    <Form.Select>

                                    </Form.Select>
                                </InputGroup>
                            </div>
                        </Form.Group>
                    </Form>
                </Modal.Body>
            </Modal>
        </>
    )
}

export default OptionModal;