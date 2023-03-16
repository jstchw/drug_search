import React from 'react';
import { Form, InputGroup, Button, OverlayTrigger, Popover } from 'react-bootstrap';
import { getDrugEventsSearch } from '../services/FDA_Request';
import SearchResultList from './SearchResultList';
import LoadingPlaceholder from './LoadingPlaceholder';
import useSearchPlaceholder from "../hooks/useSearchPlaceholder";
import { FilterLeft as FilterLeftIcon } from 'react-bootstrap-icons'
import './SearchBox.css'

const SearchBox = () => {
    const [searchTerm, setSearchTerm] = React.useState(null)

    const [searchType, setSearchType] = React.useState(null)
    const [selectedSearchTypeIndex, setSelectedSearchTypeIndex] = React.useState(0);

    const [searchResults, setSearchResults] = React.useState(null)
    const [isLoading, setIsLoading] = React.useState(false)

    const currentSearchPlaceholder = useSearchPlaceholder(3000)

    const handleSearch = async (e) => {
        e.preventDefault()
        if (searchTerm) {
            setIsLoading(true)
            setSearchResults(null)
            const data = await getDrugEventsSearch(searchTerm, searchType)
            setSearchResults(data)
            setIsLoading(false)
        }
    }

    const handleInputChange = (e) => {
        setSearchTerm(e.target.value)
    }

    const handleSearchTypeChange = (searchTypeValue, index) => {
        setSearchType(searchTypeValue)
        setSelectedSearchTypeIndex(index)
    }

    const popover = (
        <Popover id="popover-basic" className="searchbox-popover">
            <Popover.Header as="h3">Filters</Popover.Header>
            <Popover.Body>
                <p
                    className={`search-type-option${selectedSearchTypeIndex === 0 ? ' selected' : ''}`}
                    onClick={() =>
                        handleSearchTypeChange('patient.drug.openfda.brand_name', 0)
                    }
                >
                    Drug Name
                </p>
                <p
                    className={`search-type-option${selectedSearchTypeIndex === 1 ? ' selected' : ''}`}
                    onClick={() =>
                        handleSearchTypeChange('patient.reaction.reactionmeddrapt', 1)
                    }
                >
                    Side Effect
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
                            placeholder={currentSearchPlaceholder}
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

            {isLoading ? <LoadingPlaceholder/> : searchResults && <SearchResultList searchResults={searchResults} />}
        </div>
    );
}

export default SearchBox;