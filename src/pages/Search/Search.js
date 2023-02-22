import { Container, Row } from 'react-bootstrap';
import Header from '../../components/Header'
import SearchBox from '../../components/SearchBox'

const Search = () => {
    return (
        <>
            <Header />
            <Container className={'d-flex align-items-center justify-content-center h-100'}>
                <Row className={'mx-auto text-center'}>
                    <SearchBox />
                </Row>
            </Container>
        </>
    );
}

export default Search;