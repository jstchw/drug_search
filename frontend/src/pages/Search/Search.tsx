import { Row, Nav } from 'react-bootstrap';
import Header from '../../components/Header/Header';
import SearchBox from '../../components/SearchBox/SearchBox';
import { analysisFrontendUrl } from '../../constants';



const Search = () => {
  return (
    <div className="min-vh-100 d-flex flex-column justify-content-center align-items-center">
      <Row className={'mb-2'}>
        <Header />
      </Row>
      <Row className={'mt-1 mb-3'}>
        <Nav className={'p-0'} variant={'underline'}>
          <Nav.Item>
            <Nav.Link active className={'text-primary'}>Search</Nav.Link>
          </Nav.Item>
          <Nav.Item>
            <Nav.Link href={analysisFrontendUrl}>Analysis</Nav.Link>
          </Nav.Item>
        </Nav>
      </Row>
      <Row className={'mb-4 text-center'}>
        <SearchBox />
      </Row>
    </div>
  );
};

export default Search;
