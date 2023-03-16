import { Container, Row } from 'react-bootstrap';
import Header from '../../components/Header'
import './Account.css'

const Account = () => {
    return (
        <div>
            <Header />
            <Container>
                <Row>
                    <h1>Account</h1>
                </Row>
            </Container>
        </div>
    );
}

export default Account;