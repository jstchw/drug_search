import React from 'react';
import { Container, Row, Col } from 'react-bootstrap';
import Header from '../../components/Header/Header'
import SearchBox from '../../components/SearchBox/SearchBox'
import SearchResultObject from "../../components/SearchResultObject/SearchResultObject";
import LoadingPlaceholder from "../../components/SearchBox/LoadingPlaceholder";
import { getDrugEventsSearch } from "../../services/FDA_Request";
import './Search.css'

const Search = () => {
    const [searchResults, setSearchResults] = React.useState(null)
    const [isLoading, setIsLoading] = React.useState(false)
    const [searchError, setSearchError] = React.useState(false)

    const [currentSearchTerm, setCurrentSearchTerm] = React.useState('')

    const resetSearch = () => {
        setSearchResults(null)
        setIsLoading(false)
        setSearchError(false)
        setCurrentSearchTerm('')
    }

    const handleSearch = async (searchType, searchTerm) => {
        setSearchError(false)
        if(searchTerm) {
            setIsLoading(true)
            setSearchResults(null)
            const searchResults = await getDrugEventsSearch(searchTerm, searchType)
            setIsLoading(false)
            if(searchResults.result) {
                setSearchResults(searchResults)
                setCurrentSearchTerm(searchTerm)
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
                        {searchResults &&
                            <SearchResultObject
                                searchResults={searchResults}
                                searchTerm={currentSearchTerm}
                            />
                        }
                    </Col>
                </Row>
            </Container>
        </div>
    );
}

export default Search;