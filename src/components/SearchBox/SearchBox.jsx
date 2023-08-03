import React, {useEffect} from 'react';
import {Button, Form, InputGroup, Spinner} from 'react-bootstrap';
import useSearchPlaceholder from "../../hooks/useSearchPlaceholder";
import {FilterLeft as FilterLeftIcon, Search as SearchIcon} from 'react-bootstrap-icons'
import { Dropdown } from 'react-bootstrap'
import './SearchBox.css'
import OptionModal from "../OptionModal/OptionModal";
import { useSuggestions } from "../../hooks/useSuggestions";

// Options used in the Fuse.js search library
const fuseOptions = {
    keys : ['name'],
    threshold: 0.3,
}

const SearchBox = (props) => {
    const [inputValue, setInputValue] = React.useState('')

    const [showOptionModal, setShowOptionModal] = React.useState(false)

    const handleCloseOptionModal = () => setShowOptionModal(false)
    const handleShowOptionModal = () => setShowOptionModal(true)

    // Timeout for the placeholder animation
    const currentSearchPlaceholder = useSearchPlaceholder(3000, props.searchOptions.searchBy)

    // When the searched element doesn't exist in FDA's database
    const errorBox = props.searchError ? {borderColor: 'red'} : {}
    const [errorAnimation, setErrorAnimation] = React.useState(0)
    const inputRef = React.useRef(null)

    // Elements for the suggestion mechanism
    const [suggestions, setSuggestions] = React.useState([])
    const [dropdownOpen, setDropdownOpen] = React.useState(false)
    // This ref is attached to div wrapping the dropdown menu since bootstrap doesn't support refs
    const dropdownRef = React.useRef(null)

    // Handling the key control of the dropdown menu
    const [selectedSuggestionIndex, setSelectedSuggestionIndex] = React.useState(-1)

    // Clearing the input value when the additional search is closed
    useEffect(() => {
        if(!props.showAdditionalSearch && !props.isMainSearch) {
            setInputValue('')
        }
    }, [props.showAdditionalSearch, props.isMainSearch])

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
                setDropdownOpen(false);
                if (!props.showAdditionalSearch) {
                    handleSearch(null, selectedSuggestion.item)
                }
            } else {
                if (!props.showAdditionalSearch) {
                    handleSearch(null, inputValue)
                }
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

    const fuse = useSuggestions(`${window.REACT_APP_API_URL}/get_suggestions`, fuseOptions)


    // Handling clicking outside the dropdown menu
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
            props.onSearch(props.searchOptions, valueToSearch)

            if (props.searchError) {
                setErrorAnimation((prevCount) => prevCount + 1)
            }
            setDropdownOpen(false)
            // Reset the selected suggestion index
            setSelectedSuggestionIndex(-1)
        }
    }

    const { onInputChange } = props

    useEffect(() => {
        onInputChange(inputValue)
    }, [inputValue, onInputChange])

    const handleInputChange = (e) => {
        if(props.isMainSearch) props.setSearchError(false)
        setInputValue(e.target.value)

        // Suggestions should only work when a substance name or a generic name is searched for
        if(props.searchOptions.searchBy === 'patient.drug.activesubstance.activesubstancename' ||
            props.searchOptions.searchBy === 'patient.drug.openfda.generic_name') {

            if (fuse && inputValue.length >= 3 && inputValue.length <= 20) {
                const timeout = setTimeout(() => {
                    const suggestion = fuse.search(inputValue)
                    setSuggestions(suggestion.slice(0, 5))

                    if (suggestion.length > 0) {
                        setDropdownOpen(true)
                        setSelectedSuggestionIndex(-1)
                    } else {
                        setDropdownOpen(false)
                    }
                }, 300)
                return () => clearTimeout(timeout)
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
                        {props.isMainSearch &&
                            <>
                            <Button variant="outline-secondary" id="button-filters" onClick={handleShowOptionModal}>
                            <FilterLeftIcon />
                            </Button>
                            <OptionModal
                                show={showOptionModal}
                                handleClose={handleCloseOptionModal}
                                searchOptions={props.searchOptions}
                                setSearchOptions={props.setSearchOptions}
                                selectedSearchOptionIndex={props.selectedSearchOptionIndex}
                                setSelectedSearchOptionIndex={props.setSelectedSearchOptionIndex}
                            />
                            </>
                        }


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
                                                props.onInputChange(suggestion.item)
                                                setDropdownOpen(false)
                                                !props.showAdditionalSearch && handleSearch(null, suggestion.item)
                                            }}
                                            active={index === selectedSuggestionIndex}
                                        >
                                            {suggestion.item}
                                        </Dropdown.Item>
                                    ))}
                                </Dropdown.Menu>
                            </Dropdown>
                        </div>
                        {props.isMainSearch &&
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
                        </Button>}
                    </InputGroup>
                </Form.Group>
            </Form>
        </div>
    )
}

export default SearchBox