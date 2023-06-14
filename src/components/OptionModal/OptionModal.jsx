import React from 'react';
import {Modal, Form, Dropdown} from 'react-bootstrap';

// Search types that can be selected in the popover
export const searchTypes = [
    {
        value: 'patient.drug.activesubstance.activesubstancename',
        index: 0,
        label: 'Active substance'
    },
    {
        value: 'patient.reaction.reactionmeddrapt.exact',
        index: 1,
        label: 'Side Effect'
    },
    {
        value: 'patient.drug.openfda.generic_name',
        index: 2,
        label: 'Generic Name'
    },
    {
        value: 'patient.drug.openfda.brand_name',
        index: 3,
        label: 'Brand Name'
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
                    Search by:
                    <Dropdown drop={'start'}>
                        <Dropdown.Toggle variant="success" id="dropdown-basic">
                            Filters
                        </Dropdown.Toggle>

                        <Dropdown.Menu>
                            {searchTypes.map((searchType, index) => (
                                <Dropdown.Item
                                    key={index}
                                    active={props.selectedSearchTypeIndex === searchType.index}
                                    onClick={() => handleSearchTypeChange(searchType.index)}
                                >
                                    {searchType.label}
                                </Dropdown.Item>
                            ))}
                        </Dropdown.Menu>
                    </Dropdown>

                    <Form>
                        <Form.Group className="mb-3">
                            <Form.Label>Search by</Form.Label>
                            <Form.Select>
                                {/* Add your options here */}
                            </Form.Select>
                        </Form.Group>

                        <Form.Group className="mb-3">
                            <Form.Label>Age range</Form.Label>
                            <Form.Select>
                                {/* Add your options here */}
                            </Form.Select>
                        </Form.Group>

                        <Form.Group className="mb-3">
                            <Form.Label>Gender</Form.Label>
                            <Form.Select>
                                {/* Add your options here */}
                            </Form.Select>
                        </Form.Group>
                    </Form>
                </Modal.Body>
            </Modal>
        </>
    )
}

export default OptionModal;