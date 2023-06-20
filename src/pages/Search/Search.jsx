import React, {useEffect} from 'react';
import { Container, Row, Col, Button } from 'react-bootstrap';
import Header from '../../components/Header/Header'
import SearchBox from '../../components/SearchBox/SearchBox'
import SearchResultObject from "../../components/SearchResultObject/SearchResultObject";
import {getDrugEventsSearch, getEventsOverTime} from "../../services/FDA_Request";
import './Search.css'
import { LoadingContext } from "../../contexts/LoadingContext";
import { searchAgeRange, searchTypes, searchSex } from "../../components/OptionModal/OptionModal";

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
        searchBy: 0,
        sex: 0,
    })

    const [searchOptions, setSearchOptions] = React.useState({
        searchBy: searchTypes[selectedSearchOptionIndex.searchBy].value,
        sex: searchSex[selectedSearchOptionIndex.sex].value,
        age: searchAgeRange[-1]
    })

    useEffect(() => {
        if(!showAdditionalSearch) {
            setAdditionalSearchTerm(null)
        }
    }, [showAdditionalSearch])

    const handleSearch = async (searchOptions, singleSearchTerm) => {
        setEventResults(null)
        setSearchError(false)

        if(singleSearchTerm) {
            setIsLoading(true)
            const searchTerm = (additionalSearchTerm) ? [mainSearchTerm, additionalSearchTerm] : singleSearchTerm

            if (mainSearchTerm === additionalSearchTerm) {
                setSearchError(true)
                setIsLoading(false)
                return
            }
            const eventAllGroups = await getDrugEventsSearch(searchTerm, searchOptions)
            const eventsOverTime = await getEventsOverTime(searchTerm, searchOptions)

            setIsLoading(false)

            if(eventAllGroups.result) {
                setEventResults(eventAllGroups)
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
        <LoadingContext.Provider value={{isLoading}}>
            <div className="min-vh-100 d-flex flex-column">
                <Header />
                <Container fluid className={'d-flex flex-grow-1 align-items-center'}>
                    <Row className="w-100 mx-auto">
                        <Col xs="auto" className="text-center my-5 mx-auto">
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

                            {showAdditionalSearch &&
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
                            }
                            {eventResults && eventsOverTime &&
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
        </LoadingContext.Provider>
    );
}

export default Search;