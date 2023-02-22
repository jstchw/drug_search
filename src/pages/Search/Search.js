import { Container, Row } from 'react-bootstrap';
import Header from '../../components/Header'
import SearchBox from '../../components/SearchBox'
import './Search.css'

const Search = () => {
    return (
        <div>
            <Header />
            <Container className={'d-flex min-vh-100 justify-content-center align-items-center my-3'}>
                <Row>
                    <SearchBox />
                </Row>
            </Container>
        </div>
    );
}

export default Search;