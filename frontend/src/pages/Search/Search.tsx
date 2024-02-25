import { Row } from 'react-bootstrap';
import Header from '../../components/Header/Header';
import SearchBox from '../../components/SearchBox/SearchBox';
import './Search.css';

const Search = () => {
  return (
    <div className="min-vh-100 d-flex flex-column justify-content-center align-items-center">
      <Row className={'mb-4'}>
        <Header />
      </Row>
      <Row className={'mb-4 text-center'}>
        <SearchBox />
      </Row>
    </div>
  );
};

export default Search;
