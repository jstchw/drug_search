import React from 'react';
import { Form, InputGroup, Button, OverlayTrigger, Popover } from 'react-bootstrap';
import { getDrugEventsSearch } from '../services/FDA_Request';
import SearchResultList from './SearchResultList';
import LoadingPlaceholder from './LoadingPlaceholder';
import useSearchPlaceholder from "../hooks/useSearchPlaceholder";
import { FilterLeft as FilterLeftIcon, Search as SearchIcon } from 'react-bootstrap-icons'
import './SearchBox.css'

const searchTypes = [
    {
        value: 'patient.drug.openfda.brand_name',
        index: 0,
        label: 'Drug Name'
    },
    {
        value: 'patient.reaction.reactionmeddrapt',
        index: 1,
        label: 'Side Effect'
    },
]

const SearchBox = () => {
    const [searchTerm, setSearchTerm] = React.useState(null)

    const [selectedSearchTypeIndex, setSelectedSearchTypeIndex] = React.useState(0);
    const [searchType, setSearchType] = React.useState(searchTypes[selectedSearchTypeIndex].value)

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

    const handleSearchTypeChange = (index) => {
        setSearchType(searchTypes[index].value)
        setSelectedSearchTypeIndex(index)
    }

    const popover = (
        <Popover id="popover-basic" className="searchbox-popover">
            <Popover.Header as="h3">Filters</Popover.Header>
            <Popover.Body>
                {searchTypes.map((searchType) => (
                    <React.Fragment key={searchType.index}>
                        <p
                            className={`search-type-option${selectedSearchTypeIndex === searchType.index ? ' selected' : ''}`}
                            onClick={() =>
                                handleSearchTypeChange(searchType.index)
                            }
                        >
                            {searchType.label}
                        </p>
                    </React.Fragment>
                ))}
            </Popover.Body>
        </Popover>
    );


    return (
        <div>
            <Form onSubmit={handleSearch}>
                <Form.Group controlId="searchTerm">
                    <InputGroup>
                        <OverlayTrigger trigger="click" placement="left" overlay={popover} rootClose={true}>
                            <Button variant="outline-secondary" id="button-addon2">
                                <FilterLeftIcon />
                            </Button>
                        </OverlayTrigger>
                        <Form.Control
                            className="search-box"
                            type="text"
                            placeholder={currentSearchPlaceholder}
                            onChange={handleInputChange}
                        />
                        <Button variant="outline-secondary" id="button-addon2" type="submit">
                            <SearchIcon />
                        </Button>
                    </InputGroup>
                </Form.Group>
            </Form>

            {isLoading ? <LoadingPlaceholder/> : searchResults && <SearchResultList searchResults={searchResults} />}
        </div>
    );
}

export default SearchBox;