import React from 'react';
import { Container, Form, InputGroup, Button, OverlayTrigger, Popover } from 'react-bootstrap';
import { getDrugEventsSearch } from '../services/FDA_Request';
import SearchResultList from './SearchResultList';
import { FilterLeft } from 'react-bootstrap-icons'

const SearchBox = () => {
    const [searchTerm, setSearchTerm] = React.useState(null)
    const [searchType, setSearchType] = React.useState('patient.drug.openfda.brand_name')
    const [limit, setLimit] = React.useState(1)
    const [searchResults, setSearchResults] = React.useState(null)

    const handleSearch = async (e) => {
        e.preventDefault()
        setSearchResults(null)
        console.log(searchType, searchTerm, limit)
        const data = await getDrugEventsSearch(searchType, searchTerm, limit)
        setSearchResults(data)
    }

    const handleInputChange = (e) => {
        setSearchTerm(e.target.value)
    }

    const handleSearchTypeChange = (e) => {
        setSearchType(e.target.value)
    }

    const handleLimitChange = (e) => {
        setLimit(e.target.value)
    }

    const popover = (
        <Popover id="popover-basic">
            <Popover.Header as="h3">Search Types</Popover.Header>
            <Popover.Body>
                <p>Brand Name</p>
                <p>Generic Name</p>
            </Popover.Body>
        </Popover>
    );


    return (
        <div>
            <Form>
                <Form.Group controlId="searchTerm">
                    <InputGroup>
                        <Form.Control
                            className="search-box"
                            type="text"
                            placeholder={"Search for a drug..."}
                            onChange={handleInputChange}
                        />
                        <OverlayTrigger trigger="click" placement="right" overlay={popover}>
                            <Button variant="outline-secondary" id="button-addon2"><FilterLeft /></Button>
                        </OverlayTrigger>
                    </InputGroup>
                </Form.Group>



                {/*<Form.Group controlId="searchType">*/}
                {/*    <Form.Label>Search Type</Form.Label>*/}
                {/*    <Form.Control as="select" value={searchType} onChange={handleSearchTypeChange}>*/}
                {/*        <option value="patient.drug.openfda.brand_name">Brand Name</option>*/}
                {/*        <option value="patient.drug.openfda.generic_name">Generic Name</option>*/}
                {/*    </Form.Control>*/}
                {/*</Form.Group>*/}

                {/*<Form.Group controlId="limit">*/}
                {/*    <Form.Label>Limit</Form.Label>*/}
                {/*    <Form.Control type="number" min="1" max="100" value={limit} onChange={handleLimitChange}/>*/}
                {/*</Form.Group>*/}

                {/*<Container className='d-flex justify-content-center my-3'>*/}
                {/*    <button className="btn btn-primary" onClick={handleSearch}>Search</button>*/}
                {/*</Container>*/}

            </Form>

            {/*{searchResults && <SearchResultList searchResults={searchResults} />}*/}
        </div>
    );
}

export default SearchBox;