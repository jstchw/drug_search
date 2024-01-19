import React, {useEffect} from 'react';
import { Button, Form, InputGroup, Spinner } from 'react-bootstrap';
import { FilterLeft as FilterLeftIcon, Search as SearchIcon } from 'react-bootstrap-icons'
import './SearchBox.css'
import OptionModal from "../OptionModal/OptionModal";
import { useSuggestions } from "../../hooks/useSuggestions";
import { SearchOptions } from "../../types";
import {defaultSearchOptions} from "../../constants";
import Fuse, { FuseResult, Expression } from "fuse.js";
import CreatableSelect from 'react-select/creatable';
import Cookies from "js-cookie";
import { useSearch } from "../../hooks/useSearch";
//import useSearchPlaceholder from "../../hooks/useSearchPlaceholder";

type SearchBoxProps = {
    passedInput?: string[]
}


const SearchBox: React.FC<SearchBoxProps> = ({ passedInput }) => {
    const initialSearchOptionsState = () => {
        const optionsString = Cookies.get('searchOptions')
        if (optionsString) {
            try {
                return JSON.parse(optionsString)
            } catch (e) {
                console.error(e)
            }
        }
        return defaultSearchOptions
    }

    const [searchOptions, setSearchOptions] = React.useState<SearchOptions>(initialSearchOptionsState)

    useEffect(() => {
        const optionsString = Cookies.get('searchOptions')
        if (optionsString) {
            try {
                const options = JSON.parse(optionsString)
                setSearchOptions(options)
            } catch (e) {
                console.error(e)
            }
        }
    }, []);

    const [isSearching, setIsSearching] = React.useState(false)

    const [inputValue, setInputValue] = React.useState<string[]>([])
    const inputRef = React.useRef<string[]>([])

    const [showOptionModal, setShowOptionModal] = React.useState(false)

    const { executeSearch } = useSearch('SearchOptions')

    // Timeout for the placeholder animation
    //const searchPlaceholder = useSearchPlaceholder(3000, searchOptions.searchBy)

    // When the searched element doesn't exist in FDA's database
    // const errorBox = props.searchError ? {borderColor: 'red'} : {}
    // const [errorAnimation, setErrorAnimation] = React.useState(0)
    //const inputRef = React.useRef(null)

    const [fuse, setFuse] = React.useState<Fuse<string> | null>(null)
    // Elements for the suggestion mechanism
    const [suggestions, setSuggestions] = React.useState<FuseResult<string>[]>([])

    // React.useEffect(() => {
    //    if (props.searchError) {
    //        inputRef.current.classList.add('shake')
    //        const timeout = setTimeout(() => {
    //            inputRef.current.classList.remove('shake')
    //        }, 300)
    //           return () => clearTimeout(timeout)
    //    }
    // }, [errorAnimation, props.searchError])

    useSuggestions(searchOptions.searchBy, setFuse)

    const handleSearch = async (e: { preventDefault: () => void; }) => {
        setIsSearching(true)
        if (e) e.preventDefault()
        executeSearch(inputValue, searchOptions)
            // if (props.searchError) {
            //     setErrorAnimation((prevCount) => prevCount + 1)
            // }
        setIsSearching(false)
    }



    const handleInputChange = (e: string | Expression) => {
        // inputRef is used for tracking the input value in real time without state delays
        // Needed to manipulate the menu (open and close) when there are no suggestions
        inputRef.current = [e.toString()]
        if (e) {
            if (searchOptions.searchBy.index === 0 ||
                searchOptions.searchBy.index === 1) {
                if (fuse) {
                    const suggestion = fuse.search(e, {limit: 5})
                    if (suggestion.length > 0) {
                        setSuggestions(suggestion)
                    }
                }
            } else {
                    /* If the search type is side effects, it will update the input value every time the user types
                       This is needed because we don't have suggestions for side effects yet (04/01/2024)
                    */
                setInputValue([e.toString()])
            }
        } else {
            setSuggestions([])
        }
    }


    return (
        <div className={'d-flex justify-content-center'}>
            <Form onSubmit={handleSearch}>
                <Form.Group controlId="searchTerm">
                    <InputGroup>
                        <>
                        <Button variant="outline-secondary" id="button-filters" onClick={() => setShowOptionModal(true)}>
                        <FilterLeftIcon />
                        </Button>
                        <OptionModal
                            showOptionModal={showOptionModal}
                            setShowOptionModal={setShowOptionModal}
                            searchOptions={searchOptions}
                            setSearchOptions={setSearchOptions}
                        />
                        </>

                        <CreatableSelect
                            className={'creatable-select'}
                            isMulti
                            // display the passed input value if it exists
                            defaultValue={passedInput ? passedInput.map((item) => ({value: item, label: item})) : undefined}
                            menuIsOpen={suggestions.length > 0 || (searchOptions.searchBy.index === 2 && inputRef.current.some((item) => item.length > 0))}
                            onInputChange={(e) => handleInputChange(e)}
                            onChange={(e) => setInputValue(e.map((item) => item.value))}
                            options={suggestions.map((suggestion) => ({value: suggestion.item, label: suggestion.item}))}
                            components={{
                                DropdownIndicator: () => null,
                                IndicatorSeparator: () => null,
                            }}
                            styles={{
                                // Needs fixing !!!
                                control: (provided) => ({...provided, backgroundColor: 'transparent'}),
                            }}
                        />

                        <Button variant="outline-primary" id="button-submit" type="submit">
                        {isSearching ? (
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
    )
}

export default SearchBox