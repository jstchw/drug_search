import React from 'react';
import {Button, Form, InputGroup, Spinner} from 'react-bootstrap';
import useSearchPlaceholder from "../../hooks/useSearchPlaceholder";
import {FilterLeft as FilterLeftIcon, Search as SearchIcon} from 'react-bootstrap-icons'
import Fuse from 'fuse.js'
import { Dropdown } from 'react-bootstrap'
import './SearchBox.css'
import OptionModal from "../OptionModal/OptionModal";
import {searchTypes} from "../OptionModal/OptionModal";

// Options used in the Fuse.js search library
const fuseOptions = {
    keys : ['name'],
    threshold: 0.3,
}

const SearchBox = (props) => {
    const [inputValue, setInputValue] = React.useState('')

    // Search type selection
    const [selectedSearchTypeIndex, setSelectedSearchTypeIndex] = React.useState(0);
    const [searchType, setSearchType] = React.useState(searchTypes[selectedSearchTypeIndex].value)
    const [showOptionModal, setShowOptionModal] = React.useState(false)

    const handleCloseOptionModal = () => setShowOptionModal(false)
    const handleShowOptionModal = () => setShowOptionModal(true)

    // Timeout for the placeholder animation
    const currentSearchPlaceholder = useSearchPlaceholder(3000)

    // When the searched element doesn't exist in FDA's database
    const errorBox = props.searchError ? {borderColor: 'red'} : {}
    const [errorAnimation, setErrorAnimation] = React.useState(0)
    const inputRef = React.useRef(null)

    // Elements for the suggestion mechanism
    const [fuse, setFuse] = React.useState(null)
    const [suggestions, setSuggestions] = React.useState([])
    const [dropdownOpen, setDropdownOpen] = React.useState(false)
    // This ref is attached to div wrapping the dropdown menu since bootstrap doesn't support refs
    const dropdownRef = React.useRef(null)

    // Handling the key control of the dropdown menu
    const [selectedSuggestionIndex, setSelectedSuggestionIndex] = React.useState(-1)

    // The function allows to browse through suggestions using arrow keys
    const handleKeyDown = (e) => {
        if (e.keyCode === 40) {
            e.preventDefault()
            setSelectedSuggestionIndex((prevIndex) => (prevIndex + 1) % suggestions.length)
        } else if (e.keyCode === 38) {
            e.preventDefault()
            if (selectedSuggestionIndex <= 0) {
                // Set the suggestion index to the last element of the array
                setSelectedSuggestionIndex(suggestions.length - 1)
            } else {
                setSelectedSuggestionIndex((prevIndex) => (prevIndex - 1) % suggestions.length)
            }
        } else if (e.keyCode === 13) {
            e.preventDefault()
            if (suggestions.length > 0 && selectedSuggestionIndex >= 0) {
                const selectedSuggestion = suggestions[selectedSuggestionIndex];
                setInputValue(selectedSuggestion.item);
                setSelectedSuggestionIndex(-1);
                handleSearch(null, selectedSuggestion.item);
            } else {
                handleSearch(null, inputValue);
            }
        }
    }

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
                // Reset the selected suggestion index
            } catch (error) {
                console.error('Error fetching suggestions:', error)
            }
        }
        fetchSuggestions()
    }, [])


    React.useEffect(() => {
        if (dropdownOpen) {
            document.addEventListener('mousedown', handleClickOutside)
        } else {
            document.removeEventListener('mousedown', handleClickOutside)
        }
        return () => {
            document.removeEventListener('mousedown', handleClickOutside)
        }
    }, [dropdownOpen])

    const handleClickOutside = (e) => {
        if (dropdownRef.current && !dropdownRef.current.contains(e.target)) {
            setDropdownOpen(false)
        }
    }

    const handleSearch = async (e, newSearchValue) => {
        if (e) e.preventDefault()

        const valueToSearch = newSearchValue || inputValue

        if (valueToSearch) {
            props.onSearch(searchType, valueToSearch)
            if (props.searchError) {
                setErrorAnimation((prevCount) => prevCount + 1)
            }
            setDropdownOpen(false)
            // Reset the selected suggestion index
            setSelectedSuggestionIndex(-1)
        }
    }

    const handleInputChange = (e) => {
        props.setSearchError(false)
        const inputValue = e.target.value
        setInputValue(inputValue)

        // Suggestions should only work when a substance name or a generic name is searched for
        if(searchType === 'patient.drug.activesubstance.activesubstancename' ||
            searchType === 'patient.drug.openfda.generic_name') {

            if (fuse && inputValue.length >= 3) {
                const suggestion = fuse.search(inputValue)
                setSuggestions(suggestion.slice(0, 5))

                if (suggestion.length > 0) {
                    setDropdownOpen(true)
                    setSelectedSuggestionIndex(-1)
                } else {
                    setDropdownOpen(false)
                }
            } else {
                setSuggestions([])
                setDropdownOpen(false)
            }
        }
    }


    return (
        <div className={'d-flex justify-content-center'}>
            <Form onSubmit={handleSearch}>
                <Form.Group controlId="searchTerm">
                    <InputGroup>
                        <Button variant="outline-secondary" id="button-filters" onClick={handleShowOptionModal}>
                            <FilterLeftIcon />
                        </Button>
                        <OptionModal show={showOptionModal} handleClose={handleCloseOptionModal}
                                     selectedSearchTypeIndex={selectedSearchTypeIndex} setSearchType={setSearchType}
                                     setSelectedSearchTypeIndex={setSelectedSearchTypeIndex}/>

                        <div ref={dropdownRef}>
                            <Dropdown show={dropdownOpen}>
                                <Form.Control
                                    className={`search-box ${props.searchError && errorAnimation ? 'shake' : ''} rounded-0`}
                                    type="text"
                                    placeholder={currentSearchPlaceholder}
                                    onChange={handleInputChange}
                                    onKeyDown={handleKeyDown}
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
                                            active={index === selectedSuggestionIndex}
                                        >
                                            {suggestion.item}
                                        </Dropdown.Item>
                                    ))}
                                </Dropdown.Menu>
                            </Dropdown>
                        </div>
                        <Button variant="outline-primary" id="button-submit" type="submit">
                            {props.loadingSpinner ? (
                                <>
                                    <Spinner
                                        as="span"
                                        animation="border"
                                        size="sm"
                                        role="status"
                                        aria-hidden="true"
                                    />
                                </>
                            ) : (
                                <>
                                    <SearchIcon />
                                </>
                            )}
                        </Button>
                    </InputGroup>
                </Form.Group>
            </Form>
        </div>
    );
}

export default SearchBox;