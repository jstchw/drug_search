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

    React.useEffect(() => {
        console.log('searchError changed: ', searchError)
    }, [searchError])

    const handleSearch = async (searchTerm, searchType) => {
        setSearchError(false)
        if(searchTerm) {
            setIsLoading(true)
            setSearchResults(null)
            const searchResults = await getDrugEventsSearch(searchTerm, searchType)
            setIsLoading(false)
            if(searchResults.result) {
                setSearchResults(searchResults)
            } else {
                setSearchError(true)
            }
        }
    }


    return (
        <div className="min-vh-100 d-flex flex-column">
            <Header />
            <Container fluid className={'d-flex flex-grow-1 align-items-center'}>
                <Row className="w-100 justify-content-center">
                    <Col xs="auto" className="text-center">
                        <SearchBox
                            onSearch={handleSearch}
                            searchError={searchError}
                            setSearchError={setSearchError}
                        />
                        {isLoading ? (
                            <LoadingPlaceholder />
                        ) : (
                            searchResults && <SearchResultObject searchResults={searchResults} />
                        )}
                    </Col>
                </Row>
            </Container>
        </div>
    );
}

export default Search;