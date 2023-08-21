import { Container, Navbar, Badge } from 'react-bootstrap';
import { Search as SearchIcon } from 'react-bootstrap-icons';

const Header = () => {
    return (
        // Navbar containing the links to all the pages of the webapp
        <Navbar className="header-full-width">
            <Container className={'d-flex justify-content-center'}>
                <Badge className={'p-2'} style={{cursor: "pointer"}} onClick={e => window.location.href = '/'}>
                    <span className={'fs-5'}><SearchIcon/>DrugSearch</span>
                </Badge>
            </Container>
        </Navbar>
    );
}

export default Header;