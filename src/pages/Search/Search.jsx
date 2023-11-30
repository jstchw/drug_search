import React, {useEffect} from 'react'
import {Container, Row, Col, Button} from 'react-bootstrap'
import Header from '../../components/Header/Header'
import SearchBox from '../../components/SearchBox/SearchBox'
import SearchResultObject from "../../components/SearchResultObject/SearchResultObject"
import {getDrugsFromEvents, getEventsFromDrugs, getEventsOverTime} from "../../services/FDA_Request"
import './Search.css'
import {searchAgeRange, searchTypes, searchSex, searchCountry} from "../../components/OptionModal/OptionModal"
import { CSSTransition } from "react-transition-group"
import Cookies from "js-cookie"
import {ClockHistory} from "react-bootstrap-icons"
import updateSearchHistory from "../../services/updateSearchHistory"
import {SearchHistoryElement} from "../../types";


const Search = () => {
    const [eventResults, setEventResults] = React.useState(null)
    const [eventsOverTime, setEventsOverTime] = React.useState(null)


    const [isLoading, setIsLoading] = React.useState(false)
    const [searchError, setSearchError] = React.useState(false)

    const [currentSearchTerm, setCurrentSearchTerm] = React.useState('')
    const [mainSearchTerm, setMainSearchTerm] = React.useState('')
    const [additionalSearchTerm, setAdditionalSearchTerm] = React.useState('')

    const [showAdditionalSearch, setShowAdditionalSearch] = React.useState(false)

    const [selectedSearchOptionIndex, setSelectedSearchOptionIndex] = React.useState({
        searchBy: searchTypes.findIndex(({ value }) => value === (Cookies.get('searchBy') || searchTypes[0].value)),
        sex: searchSex.findIndex(({ value }) => value === (Cookies.get('lastSelectedSex') || searchSex[0].value)),
        country: searchCountry.findIndex(({ value }) => value === (Cookies.get('lastSelectedCountry') || 'US')),
    })


    const [searchOptions, setSearchOptions] = React.useState({
        searchBy: (selectedSearchOptionIndex.searchBy >= 0 && selectedSearchOptionIndex.searchBy < searchTypes.length)
            ? searchTypes[selectedSearchOptionIndex.searchBy].value : 'patient.drug.openfda.generic_name',
        sex: searchSex[selectedSearchOptionIndex.sex].value,
        age: searchAgeRange[-1],
        country: searchCountry[selectedSearchOptionIndex.country].value,
    })

    const [searchHistory, setSearchHistory] = React.useState(JSON.parse(Cookies.get('searchHistory') || '[]'))


    useEffect(() => {
        if(!showAdditionalSearch) {
            setAdditionalSearchTerm(null)
        }
    }, [showAdditionalSearch])

    // Manipulate the cookies and display of drug search terms
    useEffect(() => {
        if (currentSearchTerm) {
            let capitalizedTerms = currentSearchTerm.map((term: string) => term.charAt(0).toUpperCase() + term.slice(1));
            document.title = capitalizedTerms.join(' & ') + ' - DrugSearch';
            updateSearchHistory(capitalizedTerms, searchHistory, setSearchHistory, searchOptions)
        }
    }, [currentSearchTerm]);


    const handleSearch = async (searchOptions, singleSearchTerm) => {

        let events

        setEventResults(null)
        setSearchError(false)


        // This is a very, very bad implementation but it works
        // 100% needs to be refactored
        // When a search is submitted from a cookie button - the search term is an array from the beginning
        // When a search is submitted from the search box - the search term needs to be formed into an array and checked if other search term is present
        if(singleSearchTerm) {
            setIsLoading(true)

            let searchTerm
            if(typeof singleSearchTerm === 'string') {
                searchTerm = additionalSearchTerm ? [mainSearchTerm, additionalSearchTerm] : [singleSearchTerm]
            } else {
                searchTerm = singleSearchTerm
            }

            if (mainSearchTerm === additionalSearchTerm) {
                setSearchError(true)
                setIsLoading(false)
                return
            }

            // Search for an ADE and get drugs
            if (searchOptions.searchBy === searchTypes[2].value) {
                events = await getDrugsFromEvents(searchTerm, searchOptions)
            } else {
                events = await getEventsFromDrugs(searchTerm, searchOptions)
            }

            const eventsOverTime = await getEventsOverTime(searchTerm, searchOptions)

            setIsLoading(false)

            if(events.result) {
                setEventResults(events)
                setCurrentSearchTerm(searchTerm)
            } else {
                setSearchError(true)
            }

            if(eventsOverTime.result) {
                setEventsOverTime(eventsOverTime)
            } else {
                setSearchError(true)
            }
        }
    }

    return (
        <div className="min-vh-100 d-flex flex-column">
            <Container fluid className={'d-flex flex-grow-1 align-items-center'}>
                <Row className="w-100 mx-auto">
                    <Col xs="auto" className="text-center mb-5 mt-3 mx-auto">
                        <Row className={'mb-4'}>
                            <Container className={'mb-4 mt-3 w-auto'}>
                                <Header />
                            </Container>
                            <SearchBox
                                onSearch={handleSearch}
                                searchError={searchError}
                                setSearchError={setSearchError}
                                loadingSpinner={isLoading}
                                isMainSearch={true}
                                showAdditionalSearch={showAdditionalSearch}
                                onInputChange={setMainSearchTerm}
                                searchOptions={searchOptions}
                                setSearchOptions={setSearchOptions}
                                selectedSearchOptionIndex={selectedSearchOptionIndex}
                                setSelectedSearchOptionIndex={setSelectedSearchOptionIndex}
                            />
                        </Row>

                        {(showAdditionalSearch) ? (
                            <Button variant={'primary'} onClick={() => setShowAdditionalSearch(prevState => !prevState)}>
                                -
                            </Button>
                        ) : (
                            <Button variant={'outline-primary'}
                                    onClick={() => {
                                        setAdditionalSearchTerm(null)
                                        setShowAdditionalSearch(prevState => !prevState)
                                    }}>
                                +
                            </Button>
                        )}

                            <CSSTransition
                                in={showAdditionalSearch}
                                timeout={300}
                                classNames="fade"
                                unmountOnExit
                            >
                                <Row className={'mt-4'}>
                                    <SearchBox
                                        onSearch={handleSearch}
                                        isMainSearch={false}
                                        showAdditionalSearch={showAdditionalSearch}
                                        onInputChange={setAdditionalSearchTerm}
                                        searchOptions={searchOptions}
                                        setSearchOptions={setSearchOptions}
                                        selectedSearchOptionIndex={selectedSearchOptionIndex}
                                        setSelectedSearchOptionIndex={setSelectedSearchOptionIndex}
                                    />
                                </Row>
                            </CSSTransition>

                        {/* Recent searches powered by cookies with js-cookie */}

                        {!eventResults && searchHistory.length > 0 &&
                                <Container>
                                    <Container className={'d-flex align-items-center justify-content-center opacity-75 lead mt-4'}>
                                        <ClockHistory className={'mx-1'} />
                                        <span>Recent Searches</span>
                                    </Container>
                                    <Container className={'mt-2 d-flex flex-column align-items-center opacity-75'}>
                                        {searchHistory.map((item: SearchHistoryElement) => {
                                        return (
                                            <Button
                                                key={item.terms.join(' & ')}
                                                variant={'link'}
                                                className={'mb-2 recent-search'}
                                                onClick={() => {
                                                    handleSearch(item.options, item.terms)
                                                    setSearchOptions(item.options)
                                                }}
                                            >
                                                {item.terms.join(' & ') + ' • ' + searchTypes.find(({ value }) => value === item.options.searchBy).label}
                                            </Button>
                                        )
                                    })}</Container>
                                </Container>
                        }

                        {eventResults && eventsOverTime && currentSearchTerm && searchOptions &&
                            <SearchResultObject
                                searchResults={eventResults}
                                searchTerm={currentSearchTerm}
                                eventsOverTime={eventsOverTime}
                                showAdditionalSearch={showAdditionalSearch}
                                searchOptions={searchOptions}
                            />
                        }
                    </Col>
                </Row>
            </Container>
        </div>
    );
}

export default Search;