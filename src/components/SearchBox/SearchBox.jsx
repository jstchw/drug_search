import React from 'react';
import {Button, Form, InputGroup, OverlayTrigger, Popover} from 'react-bootstrap';
import useSearchPlaceholder from "../../hooks/useSearchPlaceholder";
import {FilterLeft as FilterLeftIcon, Search as SearchIcon} from 'react-bootstrap-icons'
import Fuse from 'fuse.js'
import { Dropdown } from 'react-bootstrap'
import './SearchBox.css'

// Search types that can be selected in the popover
const searchTypes = [
    {
        value: 'patient.drug.openfda.generic_name',
        index: 0,
        label: 'Substance Name'
    },
    {
        value: 'patient.reaction.reactionmeddrapt.exact',
        index: 1,
        label: 'Side Effect'
    },
]

// Options used in the Fuse.js search library
const fuseOptions = {
    keys : ['name'],
    threshold: 0.3,
}

const SearchBox = (props) => {
    const [inputValue, setInputValue] = React.useState('')

    const [selectedSearchTypeIndex, setSelectedSearchTypeIndex] = React.useState(0);
    const [searchType, setSearchType] = React.useState(searchTypes[selectedSearchTypeIndex].value)

    const currentSearchPlaceholder = useSearchPlaceholder(3000)

    const errorBox = props.searchError ? {borderColor: 'red'} : {}
    const [errorAnimation, setErrorAnimation] = React.useState(0)
    const inputRef = React.useRef(null)

    const [fuse, setFuse] = React.useState(null)
    const [suggestions, setSuggestions] = React.useState([])
    const [dropdownOpen, setDropdownOpen] = React.useState(false)

    React.useEffect(() => {
       if (props.searchError) {
           inputRef.current.classList.add('shake')
           const timeout = setTimeout(() => {
               inputRef.current.classList.remove('shake')
           }, 300)
              return () => clearTimeout(timeout)
       }
    }, [errorAnimation, props.searchError])

    React.useEffect(() => {
        async function fetchSuggestions(e) {
            try {
                const response = await fetch('http://localhost:16000/api/get_suggestions')
                const data = await response.json()
                setFuse(new Fuse(data, fuseOptions))
            } catch (error) {
                console.error('Error fetching suggestions:', error)
            }
        }
        fetchSuggestions()
    }, [])

    const handleSearch = async (e, newSearchValue) => {
        if (e) e.preventDefault()

        const valueToSearch = newSearchValue || inputValue

        if (valueToSearch) {
            props.onSearch(searchType, valueToSearch)
            if (props.searchError) {
                setErrorAnimation((prevCount) => prevCount + 1)
            }
            setDropdownOpen(false)
        }
    }

    const handleInputChange = (e) => {
        props.setSearchError(false)
        const inputValue = e.target.value
        setInputValue(inputValue)

        if (fuse && inputValue.length >= 3) {
            const suggestion = fuse.search(inputValue)
            console.log(suggestion)
            setSuggestions(suggestion.slice(0, 5))

            if (suggestion.length > 0) {
                setDropdownOpen(true)
            } else {
                setDropdownOpen(false)
            }
        } else {
            setSuggestions([])
            setDropdownOpen(false)
        }
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
                            <Button variant="outline-secondary" id="button-filters">
                                <FilterLeftIcon />
                            </Button>
                        </OverlayTrigger>
                        <Dropdown show={dropdownOpen}>
                            <Form.Control
                                className={`search-box ${props.searchError && errorAnimation ? 'shake' : ''}`}
                                type="text"
                                placeholder={currentSearchPlaceholder}
                                onChange={handleInputChange}
                                style={errorBox}
                                ref={inputRef}
                                value={inputValue}
                                autoComplete={'off'}
                                autoCorrect={'off'}
                                autoCapitalize={'off'}
                                spellCheck={'false'}
                            />
                            <Dropdown.Menu className="search-suggestions">
                                {suggestions.map((suggestion, index) => (
                                    <Dropdown.Item
                                        key={index}
                                        onClick={() => {
                                            setInputValue(suggestion.item)
                                            handleSearch(null, suggestion.item)
                                        }}
                                    >
                                        {suggestion.item}
                                    </Dropdown.Item>
                                ))}
                            </Dropdown.Menu>
                        </Dropdown>
                        <Button variant="outline-primary" id="button-submit" type="submit">
                            <SearchIcon />
                        </Button>
                    </InputGroup>
                </Form.Group>
            </Form>
        </div>
    );
}

export default SearchBox;