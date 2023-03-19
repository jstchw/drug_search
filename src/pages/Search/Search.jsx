import { Container, Row, Col } from 'react-bootstrap';
import Header from '../../components/Header/Header'
import SearchBox from '../../components/SearchBox/SearchBox'
import './Search.css'

const Search = () => {
    return (
        <div className="min-vh-100 d-flex flex-column">
            <Header />
            <Container fluid className={'d-flex flex-grow-1 align-items-center'}>
                <Row className="w-100 justify-content-center">
                    <Col xs="auto" className="text-center">
                        <SearchBox />
                    </Col>
                </Row>
            </Container>
        </div>
    );
}

export default Search;