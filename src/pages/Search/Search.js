import Container from 'react-bootstrap/Container';
import Header from '../../components/Header'
import SearchBox from '../../components/SearchBox'

const Search = () => {
    return (
        <>
            <Header />
            <Container>
                <h1>Search</h1>
                <SearchBox />
            </Container>
        </>
    );
}

export default Search;