import React, {useEffect} from 'react'
import { Container, Row, Col, Button } from 'react-bootstrap'
import Header from '../../components/Header/Header'
import SearchBox from '../../components/SearchBox/SearchBox'
import SearchResultObject from "../../components/SearchResultObject/SearchResultObject"
import {getDrugsFromEvents, getEventsFromDrugs, getEventsOverTime} from "../../services/FDA_Request"
import './Search.css'
import {searchAgeRange, searchTypes, searchSex, searchCountry} from "../../components/OptionModal/OptionModal"
import { CSSTransition } from "react-transition-group"
import Cookies from "js-cookie"


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
        searchBy: searchTypes[selectedSearchOptionIndex.searchBy].value,
        sex: searchSex[selectedSearchOptionIndex.sex].value,
        age: searchAgeRange[-1],
        country: searchCountry[selectedSearchOptionIndex.country].value,
    })


    useEffect(() => {
        if(!showAdditionalSearch) {
            setAdditionalSearchTerm(null)
        }
    }, [showAdditionalSearch])

    const handleSearch = async (searchOptions, singleSearchTerm) => {

        let events

        setEventResults(null)
        setSearchError(false)

        if(singleSearchTerm) {
            setIsLoading(true)
            const searchTerm = (additionalSearchTerm) ? [mainSearchTerm, additionalSearchTerm] : [singleSearchTerm]


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
                            <Header />
                        </Row>
                        <Row className={'mb-4'}>
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