import React from 'react';
import {Modal, Form, Badge} from 'react-bootstrap';

// Search types that can be selected in the popover
export const searchTypes = [
    {
        value: 'patient.drug.activesubstance.activesubstancename',
        index: 0,
        label: 'Active substance'
    },
    {
        value: 'patient.drug.openfda.brand_name',
        index: 1,
        label: 'Brand Name'
    },
    {
        value: 'patient.reaction.reactionmeddrapt.exact',
        index: 2,
        label: 'Side Effect'
    },
]

export const searchSex = [
    {
        value: 'all',
        index: 0,
        label: 'All'
    },
    {
        value: 'male',
        index: 1,
        label: 'Male'
    },
    {
        value: 'female',
        index: 2,
        label: 'Female'
    }
]



const OptionModal = (props) => {

    const handleSearchTypeChange = (index) => {
        props.setSearchType(searchTypes[index].value)
        props.setSelectedSearchTypeIndex(index)
    }

    return (
        <>
            <Modal centered show={props.show} onHide={props.handleClose}>
                <Modal.Header closeButton>
                    <Modal.Title>Options</Modal.Title>
                </Modal.Header>
                <Modal.Body>
                    <Form>
                        <Form.Group className="mb-3">
                            <Form.Label><Badge>Search by</Badge></Form.Label>
                            <Form.Select
                                onChange={(e) => handleSearchTypeChange(e.target.value)}
                                value={props.selectedSearchTypeIndex}
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
                        </Form.Group>
                        <Form.Group className="mb-3">
                            <Form.Label><Badge>Sex</Badge></Form.Label>
                            <Form.Select>
                                {searchSex.map((searchType, index) => (
                                    <option
                                        key={index}
                                        value={searchType.index}
                                    >
                                        {searchType.label}
                                    </option>
                                ))}
                            </Form.Select>
                        </Form.Group>
                    </Form>
                </Modal.Body>
            </Modal>
        </>
    )
}

export default OptionModal;