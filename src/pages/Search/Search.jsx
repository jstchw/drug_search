import React from 'react';
import { Container, Row, Col } from 'react-bootstrap';
import Header from '../../components/Header/Header'
import SearchBox from '../../components/SearchBox/SearchBox'
import SearchResultObject from "../../components/SearchResultObject/SearchResultObject";
import {getDrugEventsSearch, getEventsOverTime} from "../../services/FDA_Request";
import './Search.css'

const Search = () => {
    const [eventResults, setEventResults] = React.useState(null)
    const [eventsOverTime, setEventsOverTime] = React.useState(null)


    const [isLoading, setIsLoading] = React.useState(false)
    const [searchError, setSearchError] = React.useState(false)

    const [currentSearchTerm, setCurrentSearchTerm] = React.useState('')

    const resetSearch = () => {
        setEventResults(null)
        setIsLoading(false)
        setSearchError(false)
        setCurrentSearchTerm('')
    }

    const handleSearch = async (searchType, searchTerm) => {
        setSearchError(false)
        if(searchTerm) {
            setIsLoading(true)
            const eventAllGroups = await getDrugEventsSearch(searchTerm, searchType)
            const eventsOverTime = await getEventsOverTime(searchTerm, searchType)
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
        <div className="min-vh-100 d-flex flex-column">
            <Header onLogoClick={resetSearch}/>
            <Container fluid className={'d-flex flex-grow-1 align-items-center'}>
                <Row className="w-100 justify-content-center">
                    <Col xs="auto" className="text-center my-5">
                        <SearchBox
                            onSearch={handleSearch}
                            searchError={searchError}
                            setSearchError={setSearchError}
                            loadingSpinner={isLoading}
                        />
                        {eventResults && eventsOverTime &&
                            <SearchResultObject
                                searchResults={eventResults}
                                searchTerm={currentSearchTerm}
                                eventsOverTime={eventsOverTime}
                            />
                        }
                    </Col>
                </Row>
            </Container>
        </div>
    );
}

export default Search;