import React from 'react';
import {Button, Form, InputGroup, OverlayTrigger, Popover} from 'react-bootstrap';
import useSearchPlaceholder from "../../hooks/useSearchPlaceholder";
import {FilterLeft as FilterLeftIcon, Search as SearchIcon} from 'react-bootstrap-icons'
import './SearchBox.css'

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

const SearchBox = (props) => {
    const [inputValue, setInputValue] = React.useState('')

    const [selectedSearchTypeIndex, setSelectedSearchTypeIndex] = React.useState(0);
    const [searchType, setSearchType] = React.useState(searchTypes[selectedSearchTypeIndex].value)

    const currentSearchPlaceholder = useSearchPlaceholder(3000)

    const errorBox = props.searchError ? {borderColor: 'red'} : {}
    const [errorAnimation, setErrorAnimation] = React.useState(0)
    const inputRef = React.useRef(null)

    React.useEffect(() => {
       if (props.searchError) {
           inputRef.current.classList.add('shake')
           const timeout = setTimeout(() => {
               inputRef.current.classList.remove('shake')
           }, 300)
              return () => clearTimeout(timeout)
       }
    }, [errorAnimation, props.searchError])

    const handleSearch = async (e) => {
        e.preventDefault()
        if (inputValue) {
            props.onSearch(searchType, inputValue)
            if (props.searchError) {
                setErrorAnimation((prevCount) => prevCount + 1)
            }
        }
    }

    const handleInputChange = (e) => {
        props.setSearchError(false)
        setInputValue(e.target.value)
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
                        <Form.Control
                            className={`search-box ${props.searchError && errorAnimation ? 'shake' : ''}`}
                            type="text"
                            placeholder={currentSearchPlaceholder}
                            onChange={handleInputChange}
                            style={errorBox}
                            ref={inputRef}
                            value={inputValue}
                        />
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