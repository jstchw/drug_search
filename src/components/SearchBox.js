import React from 'react';
import { Container, Form, InputGroup, Button, OverlayTrigger, Popover } from 'react-bootstrap';
import { getDrugEventsSearch } from '../services/FDA_Request';
import SearchResultList from './SearchResultList';
import { FilterLeft as FilterLeftIcon } from 'react-bootstrap-icons'
import './SearchBox.css'

const SearchBox = () => {
    const [searchTerm, setSearchTerm] = React.useState(null)

    const [searchType, setSearchType] = React.useState(null)
    const [selectedSearchTypeIndex, setSelectedSearchTypeIndex] = React.useState(0);

    const [limit, setLimit] = React.useState(1)
    const [searchResults, setSearchResults] = React.useState(null)

    const handleSearch = async (e) => {
        e.preventDefault()
        if (searchTerm) {
            const data = await getDrugEventsSearch(searchTerm, searchType, limit)
            setSearchResults(data)
        }
    }

    const handleInputChange = (e) => {
        setSearchTerm(e.target.value)
    }

    const handleSearchTypeChange = (searchTypeValue, index) => {
        setSearchType(searchTypeValue)
        setSelectedSearchTypeIndex(index)
    }

    const handleLimitChange = (e) => {
        setLimit(e.target.value)
    }

    const popover = (
        <Popover id="popover-basic" className="searchbox-popover">
            <Popover.Header as="h3">Search Types</Popover.Header>
            <Popover.Body>
                <p
                    className="search-type-option"
                >
                    Brand Name
                </p>
                <p
                    className={`search-type-option${selectedSearchTypeIndex === 1 ? ' selected' : ''}`}
                    onClick={() =>
                        handleSearchTypeChange('patient.drug.openfda.brand_name', 1)
                    }
                >
                    Generic Name
                </p>
            </Popover.Body>
        </Popover>
    );


    return (
        <div>
            <Form onSubmit={handleSearch}>
                <Form.Group controlId="searchTerm">
                    <InputGroup>
                        <Form.Control
                            className="search-box"
                            type="text"
                            placeholder={"Search for a drug..."}
                            onChange={handleInputChange}
                        />
                        <OverlayTrigger trigger="click" placement="right" overlay={popover} rootClose={true}>
                            <Button variant="outline-secondary" id="button-addon2">
                                <FilterLeftIcon />
                            </Button>
                        </OverlayTrigger>
                    </InputGroup>
                </Form.Group>
                <Button style={{ display: "none" }} type="submit" />
            </Form>

            {searchResults && <SearchResultList searchResults={searchResults} />}
        </div>
    );
}

export default SearchBox;